#!/usr/bin/env python3
"""Export all CardioChrom folds into a portable, deployment-only bundle.

This script must run on UA Puma beside the frozen analysis assets. It creates a
new directory and never modifies the authoritative source artifacts.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np


PROJECT = Path("/xdisk/yitan/hongbiaohuang/HFpEF/Zenodo15232790")
DEPLOYMENT_SCRIPT = Path(
    "/home/u10/hongbiaohuang/Script/fnih_fig5f_scp3342_virtual_epigenome_deployment_v1.py"
)
CELL_TYPES = (
    "vCM", "aCM", "Adipocyte", "Fibroblast", "Endothelial", "Endocardial",
    "Epicardial", "Pericyte", "Myeloid", "SM", "Lymphoid", "Neuronal",
)
MODALITIES = ("ATAC", "H3K27ac", "H3K27me3")
N_GENES = 36_100
N_ATAC = 285_873
N_HISTONE = 191_678
LATENT = 50
KNN_K = 25


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export portable CardioChrom model bundles on Puma")
    parser.add_argument("--project-root", type=Path, default=PROJECT)
    parser.add_argument("--deployment-script", type=Path, default=DEPLOYMENT_SCRIPT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/xdisk/yitan/hongbiaohuang/HFpEF/CardioChrom_model_bundle_v1"),
    )
    parser.add_argument(
        "--atac-h5ad",
        type=Path,
        help="FNIH ATAC AnnData defining the exact 285,873-feature order; auto-detected if omitted",
    )
    parser.add_argument("--skip-checksums", action="store_true")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare-common", action="store_true")
    mode.add_argument("--fold-index", type=int, choices=range(12))
    mode.add_argument("--finalize", action="store_true")
    mode.add_argument("--all", action="store_true", help="Sequential export; Slurm array mode is preferred")
    return parser.parse_args()


def import_module(path: Path):
    spec = importlib.util.spec_from_file_location("cardiochrom_frozen_deployment", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import deployment script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require_empty_output(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise SystemExit(f"Output directory is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def copy_required(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source, destination)


def feature_ids(path: Path, expected_features: int) -> list[str]:
    data = ad.read_h5ad(path, backed="r")
    try:
        if data.n_vars != expected_features:
            raise ValueError(f"{path} has {data.n_vars} features, expected {expected_features}")
        return data.var_names.astype(str).tolist()
    finally:
        data.file.close()


def find_atac_h5ad(project_root: Path) -> Path:
    matches: list[Path] = []
    for path in sorted(project_root.glob("*.h5ad")):
        if "atac" not in path.name.lower():
            continue
        try:
            data = ad.read_h5ad(path, backed="r")
            n_vars = data.n_vars
            data.file.close()
        except Exception:
            continue
        if n_vars == N_ATAC:
            matches.append(path)
    if len(matches) != 1:
        raise RuntimeError(
            "Could not uniquely auto-detect the FNIH 285,873-feature ATAC AnnData. "
            f"Candidates={matches}. Supply --atac-h5ad explicitly."
        )
    return matches[0]


def write_features(path: Path, values: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("feature_index\tfeature_id\n")
        for index, value in enumerate(values):
            handle.write(f"{index}\t{value}\n")


def array_shape(path: Path) -> list[int]:
    return list(np.load(path, mmap_mode="r", allow_pickle=False).shape)


def export_existing_bundle(
    destination: Path,
    cell_type: str,
    modality: str,
    sources: dict[str, Path],
) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    for target_name, source in sources.items():
        copy_required(source, destination / f"{target_name}.npy")
    metadata = {
        "bundle_version": "1.0.0-rc2",
        "cell_type": cell_type,
        "fold": CELL_TYPES.index(cell_type),
        "modality": modality,
        "latent_dim": LATENT,
        "knn_k": KNN_K,
        "knn_weights": "uniform",
        "knn_metric": "euclidean",
        "rna_normalization": "library-size 10000 then log1p",
        "decoder": (
            "continuous latent @ components"
            if modality == "ATAC"
            else "latent @ components; clip 0..20; expm1; normalize 10000; log1p"
        ),
        "arrays": {
            name: array_shape(destination / f"{name}.npy") for name in sources
        },
        "source_artifact_names": {name: source.name for name, source in sources.items()},
    }
    write_json(destination / "metadata.json", metadata)
    return metadata


def export_h3me3_bundle(destination: Path, cell_type: str, module) -> dict:
    destination.mkdir(parents=True, exist_ok=False)
    fold = CELL_TYPES.index(cell_type)
    route = {"fold": fold, "held_out_cell_type": cell_type, "route_name": cell_type}
    rna_components, knn, _ = module.fit_h3me3(route)
    if not hasattr(knn, "_fit_X") or not hasattr(knn, "_y"):
        raise RuntimeError("Unexpected scikit-learn KNN internals; cannot export H3K27me3 latents")
    source_dir = module.ME3_ROOT / f"fold_{fold:02d}_{cell_type}"
    np.save(destination / "rna_components.npy", np.asarray(rna_components, dtype=np.float32), allow_pickle=False)
    np.save(destination / "train_rna_latent.npy", np.asarray(knn._fit_X, dtype=np.float32), allow_pickle=False)
    np.save(destination / "train_target_latent.npy", np.asarray(knn._y, dtype=np.float32), allow_pickle=False)
    copy_required(
        source_dir / "h3k27me3_svd_components.float32.npy",
        destination / "target_components.npy",
    )
    names = ("rna_components", "train_rna_latent", "train_target_latent", "target_components")
    metadata = {
        "bundle_version": "1.0.0-rc2",
        "cell_type": cell_type,
        "fold": fold,
        "modality": "H3K27me3",
        "latent_dim": LATENT,
        "knn_k": KNN_K,
        "knn_weights": "uniform",
        "knn_metric": "euclidean",
        "rna_normalization": "library-size 10000 then log1p",
        "decoder": "latent @ components; clip 0..20; expm1; normalize 10000; log1p",
        "portable_freeze": "train-only RNA and H3K27me3 latents materialized during export",
        "arrays": {name: array_shape(destination / f"{name}.npy") for name in names},
        "source_artifact_names": {
            "rna_components": "rna_svd_components.float32.npy",
            "target_components": "h3k27me3_svd_components.float32.npy",
        },
    }
    write_json(destination / "metadata.json", metadata)
    return metadata


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_checksums(output_dir: Path) -> None:
    files = sorted(path for path in output_dir.rglob("*") if path.is_file() and path.name != "checksums.sha256")
    with (output_dir / "checksums.sha256").open("w", encoding="utf-8") as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.relative_to(output_dir)}\n")


def prepare_common(args: argparse.Namespace, project_root: Path, output_dir: Path) -> None:
    require_empty_output(output_dir)
    common = output_dir / "common"
    common.mkdir()
    genes_source = project_root / "FNIH_H3K27ac_RNA_counts_67453x36100_genes.txt"
    genes = genes_source.read_text(encoding="utf-8").splitlines()
    if len(genes) != N_GENES or len(set(genes)) != N_GENES:
        raise ValueError("Frozen gene file must contain 36,100 unique genes")
    copy_required(genes_source, common / "genes.txt")

    atac_h5ad = args.atac_h5ad.resolve() if args.atac_h5ad else find_atac_h5ad(project_root)
    h3ac_h5ad = project_root / "FNIH_Heart_pool.H3K27ac.DNA.h5ad"
    h3me3_h5ad = project_root / "FNIH_Heart_pool.H3K27me3.DNA.h5ad"
    write_features(common / "atac_features.tsv.gz", feature_ids(atac_h5ad, N_ATAC))
    h3ac_features = feature_ids(h3ac_h5ad, N_HISTONE)
    h3me3_features = feature_ids(h3me3_h5ad, N_HISTONE)
    if h3ac_features != h3me3_features:
        raise ValueError("H3K27ac and H3K27me3 feature orders differ")
    write_features(common / "histone_bins.tsv.gz", h3ac_features)
    write_json(
        output_dir / "EXPORT_IN_PROGRESS.json",
        {
            "status": "common_prepared",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "atac_feature_source_name": atac_h5ad.name,
        },
    )
    print(f"Prepared common interface: {common}")


def export_fold(project_root: Path, output_dir: Path, deployment_script: Path, fold: int) -> None:
    if not (output_dir / "common" / "genes.txt").is_file():
        raise RuntimeError("Run --prepare-common before exporting folds")
    cell_type = CELL_TYPES[fold]
    fold_dir = output_dir / "folds" / f"{fold:02d}_{cell_type}"
    if fold_dir.exists():
        raise RuntimeError(f"Fold output already exists; refusing overwrite: {fold_dir}")
    module = import_module(deployment_script)
    atac_source = project_root / "FNIH_RNA2ATAC_V2_LS_GUIDED_protocol" / f"fold_{fold:02d}_{cell_type}"
    export_existing_bundle(
        fold_dir / "ATAC",
        cell_type,
        "ATAC",
        {
            "rna_components": atac_source / "rna_svd50_components.float32.npy",
            "train_rna_latent": atac_source / "rna_svd50_train_cells.float32.npy",
            "train_target_latent": atac_source / "atac_svd50_train_cells.float32.npy",
            "target_components": atac_source / "atac_svd50_components.float32.npy",
        },
    )
    h3ac_source = project_root / "FNIH_RNA2H3K27ac_LS_representation_selection_v1" / f"fold_{fold:02d}_{cell_type}"
    export_existing_bundle(
        fold_dir / "H3K27ac",
        cell_type,
        "H3K27ac",
        {
            "rna_components": h3ac_source / "rna_svd_components.float32.npy",
            "train_rna_latent": h3ac_source / "train_rna_latent.float32.npy",
            "train_target_latent": h3ac_source / "train_target_A_latent.float32.npy",
            "target_components": h3ac_source / "target_A_svd_components.float32.npy",
        },
    )
    export_h3me3_bundle(fold_dir / "H3K27me3", cell_type, module)
    print(f"Exported fold {fold:02d} {cell_type}", flush=True)


def finalize(output_dir: Path, skip_checksums: bool) -> None:
    missing: list[str] = []
    for fold, cell_type in enumerate(CELL_TYPES):
        for modality in MODALITIES:
            path = output_dir / "folds" / f"{fold:02d}_{cell_type}" / modality / "metadata.json"
            if not path.is_file():
                missing.append(str(path))
    if missing:
        raise RuntimeError("Cannot finalize; missing bundle files:\n" + "\n".join(missing))
    manifest = {
        "name": "CardioChrom portable frozen model bundle",
        "bundle_version": "1.0.0-rc2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "cell_types": list(CELL_TYPES),
        "modalities": list(MODALITIES),
        "n_genes": N_GENES,
        "n_atac_features": N_ATAC,
        "n_histone_features": N_HISTONE,
        "latent_dim": LATENT,
        "knn_k": KNN_K,
        "n_bundles": len(CELL_TYPES) * len(MODALITIES),
        "validation_status": "pending SCP3342 frozen-output comparison",
    }
    write_json(output_dir / "bundle_manifest.json", manifest)
    in_progress = output_dir / "EXPORT_IN_PROGRESS.json"
    if in_progress.exists():
        in_progress.unlink()
    if not skip_checksums:
        write_checksums(output_dir)
    print(f"Portable bundle finalized: {output_dir}")


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    deployment_script = args.deployment_script.resolve()
    if not project_root.is_dir() or not deployment_script.is_file():
        raise FileNotFoundError("Project root or deployment script is missing")
    if not output_dir.parent.is_dir():
        raise FileNotFoundError(output_dir.parent)
    if shutil.disk_usage(output_dir.parent).free < 5 * 1024**3:
        raise RuntimeError("At least 5 GiB of free space is required for the portable bundle")
    if args.prepare_common:
        prepare_common(args, project_root, output_dir)
    elif args.fold_index is not None:
        export_fold(project_root, output_dir, deployment_script, args.fold_index)
    elif args.finalize:
        finalize(output_dir, args.skip_checksums)
    elif args.all:
        prepare_common(args, project_root, output_dir)
        for fold in range(12):
            export_fold(project_root, output_dir, deployment_script, fold)
        finalize(output_dir, args.skip_checksums)
    return 0


if __name__ == "__main__":
    sys.exit(main())
