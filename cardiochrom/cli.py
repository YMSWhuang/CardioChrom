from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .bundle import BundleRegistry
from .constants import CANONICAL_CELL_TYPES, MODALITIES
from .core import predict_latent


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="cardiochrom")
    commands = root.add_subparsers(dest="command", required=True)
    predict = commands.add_parser("predict", help="Run frozen CardioChrom inference")
    predict.add_argument("--input", type=Path, required=True)
    predict.add_argument("--model-dir", type=Path, required=True)
    predict.add_argument("--cell-type-key", required=True)
    predict.add_argument("--routing-map", type=Path)
    predict.add_argument("--layer", default=None, help="Counts layer; default uses adata.X")
    predict.add_argument("--modalities", default=",".join(MODALITIES))
    predict.add_argument("--output", type=Path, required=True)
    predict.add_argument("--min-gene-coverage", type=float, default=1.0)
    predict.add_argument("--n-jobs", type=int, default=1)
    return root


def load_routes(path: Path | None) -> dict[str, str]:
    if path is None:
        return {cell_type: cell_type for cell_type in CANONICAL_CELL_TYPES}
    routes = pd.read_csv(path, sep="\t", dtype=str)
    required = {"input_label", "canonical_cell_type"}
    if not required.issubset(routes.columns):
        raise ValueError(f"Routing map must contain columns: {sorted(required)}")
    if routes["input_label"].duplicated().any():
        raise ValueError("Routing map contains duplicate input labels")
    invalid = sorted(set(routes["canonical_cell_type"]) - set(CANONICAL_CELL_TYPES))
    if invalid:
        raise ValueError(f"Routing map contains unsupported canonical cell types: {invalid}")
    return dict(zip(routes["input_label"], routes["canonical_cell_type"]))


def align_genes(
    matrix: sp.spmatrix | np.ndarray,
    input_genes: list[str],
    model_genes: tuple[str, ...],
) -> tuple[sp.csr_matrix, float, int]:
    if len(input_genes) != len(set(input_genes)):
        raise ValueError("Input gene names must be unique")
    lookup = {gene: index for index, gene in enumerate(input_genes)}
    model_positions = np.asarray(
        [index for index, gene in enumerate(model_genes) if gene in lookup], dtype=np.int64
    )
    input_positions = np.asarray(
        [lookup[model_genes[index]] for index in model_positions], dtype=np.int64
    )
    coverage = len(model_positions) / len(model_genes)
    selected = matrix[:, input_positions]
    selected = selected.tocoo() if sp.issparse(selected) else sp.coo_matrix(selected)
    aligned = sp.coo_matrix(
        (selected.data, (selected.row, model_positions[selected.col])),
        shape=(matrix.shape[0], len(model_genes)),
        dtype=np.float32,
    ).tocsr()
    aligned.sum_duplicates()
    aligned.sort_indices()
    return aligned, coverage, len(model_genes) - len(model_positions)


def validate_raw_counts(matrix: sp.spmatrix | np.ndarray) -> None:
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size == 0:
        return
    if not np.isfinite(values).all() or np.min(values) < 0:
        raise ValueError("RNA input must contain finite, non-negative raw counts")
    if np.max(np.abs(values - np.rint(values))) > 1e-4:
        raise ValueError(
            "RNA input does not appear to contain raw counts. Supply an integer-like counts layer with --layer."
        )


def run_predict(args: argparse.Namespace) -> int:
    import anndata as ad

    registry = BundleRegistry(args.model_dir)
    data = ad.read_h5ad(args.input)
    if args.cell_type_key not in data.obs:
        raise KeyError(f"Missing obs column: {args.cell_type_key}")
    matrix = data.X if args.layer is None else data.layers[args.layer]
    matrix, coverage, missing = align_genes(matrix, data.var_names.astype(str).tolist(), registry.genes)
    validate_raw_counts(matrix)
    if coverage < args.min_gene_coverage:
        raise ValueError(
            f"Frozen-gene coverage {coverage:.6f} is below --min-gene-coverage "
            f"{args.min_gene_coverage:.6f}; missing genes={missing}"
        )

    requested = tuple(item.strip() for item in args.modalities.split(",") if item.strip())
    invalid = sorted(set(requested) - set(MODALITIES))
    if invalid:
        raise ValueError(f"Unsupported modalities: {invalid}")

    routes = load_routes(args.routing_map)
    labels = data.obs[args.cell_type_key].astype(str)
    canonical = labels.map(routes)
    args.output.mkdir(parents=True, exist_ok=True)

    outputs: dict[str, np.ndarray] = {}
    for modality in requested:
        outputs[f"{modality}_latent"] = np.full((data.n_obs, 50), np.nan, dtype=np.float32)
        outputs[f"{modality}_knn25_mean_distance"] = np.full(data.n_obs, np.nan, dtype=np.float32)

    route_counts: dict[str, int] = {}
    for cell_type in CANONICAL_CELL_TYPES:
        rows = np.flatnonzero(canonical.to_numpy(dtype=object) == cell_type)
        if not len(rows):
            continue
        route_counts[cell_type] = int(len(rows))
        counts = matrix[rows, :]
        for modality in requested:
            latent, distance = predict_latent(
                counts,
                registry.bundle(cell_type, modality),
                n_jobs=args.n_jobs,
            )
            outputs[f"{modality}_latent"][rows] = latent
            outputs[f"{modality}_knn25_mean_distance"][rows] = distance

    np.savez_compressed(args.output / "CardioChrom_latent_predictions.npz", **outputs)
    cell_table = pd.DataFrame(
        {
            "cell_id": data.obs_names.astype(str),
            "input_cell_type": labels.to_numpy(),
            "canonical_cell_type": canonical.fillna("").to_numpy(),
            "routed": canonical.notna().to_numpy(),
        }
    )
    cell_table.to_csv(args.output / "CardioChrom_cells.tsv.gz", sep="\t", index=False)
    manifest = {
        "status": "completed",
        "input": str(args.input.resolve()),
        "model_dir": str(args.model_dir.resolve()),
        "cell_type_key": args.cell_type_key,
        "layer": args.layer,
        "modalities": list(requested),
        "n_cells": int(data.n_obs),
        "n_routed_cells": int(canonical.notna().sum()),
        "route_counts": route_counts,
        "frozen_gene_coverage": coverage,
        "missing_frozen_genes": missing,
        "normalization": "library-size 10000 then log1p",
        "translator": "uniform Euclidean KNN25 in frozen 50-D latent space",
    }
    (args.output / "CardioChrom_run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return 0


def main() -> int:
    args = parser().parse_args()
    if args.command == "predict":
        return run_predict(args)
    raise RuntimeError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
