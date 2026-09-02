#!/usr/bin/env python3
"""Create a tiny synthetic CardioChrom input file.

This file is an interface smoke test only. The generated counts are random and
must not be interpreted biologically.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a tiny synthetic raw-count AnnData for CardioChrom."
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        required=True,
        help="Extracted CardioChrom model-bundle directory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("example_input.h5ad"),
        help="Output .h5ad path (default: example_input.h5ad).",
    )
    parser.add_argument("--seed", type=int, default=20260902)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    genes_path = args.model_dir / "common" / "genes.txt"
    if not genes_path.is_file():
        raise FileNotFoundError(f"Frozen gene list not found: {genes_path}")

    genes = [line.strip() for line in genes_path.read_text(encoding="utf-8").splitlines()]
    genes = [gene for gene in genes if gene]
    if not genes:
        raise ValueError(f"Frozen gene list is empty: {genes_path}")
    if len(genes) != len(set(genes)):
        raise ValueError(f"Frozen gene list contains duplicate names: {genes_path}")

    cell_types = ["vCM", "Fibroblast", "Endothelial"]
    rng = np.random.default_rng(args.seed)
    n_nonzero = min(64, len(genes))

    rows: list[int] = []
    cols: list[int] = []
    values: list[int] = []
    for row in range(len(cell_types)):
        selected = rng.choice(len(genes), size=n_nonzero, replace=False)
        rows.extend([row] * n_nonzero)
        cols.extend(selected.tolist())
        values.extend(rng.integers(1, 9, size=n_nonzero).tolist())

    counts = sparse.csr_matrix(
        (np.asarray(values, dtype=np.float32), (rows, cols)),
        shape=(len(cell_types), len(genes)),
        dtype=np.float32,
    )
    obs = pd.DataFrame(
        {"cell_type": cell_types},
        index=[f"synthetic_cell_{index + 1}" for index in range(len(cell_types))],
    )
    var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.uns["note"] = (
        "Synthetic CardioChrom interface smoke test; not biologically meaningful."
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(
        f"Wrote {args.output} with {adata.n_obs} cells, {adata.n_vars} genes, "
        f"and {counts.nnz} non-zero raw counts."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
