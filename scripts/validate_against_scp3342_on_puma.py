#!/usr/bin/env python3
"""Compare portable bundles with the authoritative frozen SCP3342 deployment."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_REPO = Path("/xdisk/yitan/hongbiaohuang/CardioChrom_release_candidate")
DEFAULT_MODEL = Path("/xdisk/yitan/hongbiaohuang/HFpEF/CardioChrom_model_bundle_v1")
DEFAULT_DEPLOYMENT_SCRIPT = Path(
    "/home/u10/hongbiaohuang/Script/fnih_fig5f_scp3342_virtual_epigenome_deployment_v1.py"
)
DEFAULT_OUTPUT = Path(
    "/xdisk/yitan/hongbiaohuang/HFpEF/CardioChrom_model_bundle_v1_SCP3342_validation.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate CardioChrom bundles against frozen SCP3342 outputs")
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--deployment-script", type=Path, default=DEFAULT_DEPLOYMENT_SCRIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sample-cells", type=int, default=32)
    parser.add_argument("--tolerance", type=float, default=2e-5)
    parser.add_argument("--n-jobs", type=int, default=4)
    return parser.parse_args()


def import_module(path: Path):
    spec = importlib.util.spec_from_file_location("cardiochrom_frozen_deployment", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import deployment script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    args = parse_args()
    sys.path.insert(0, str(args.repo_root.resolve()))
    from cardiochrom.bundle import BundleRegistry
    from cardiochrom.core import predict_latent

    module = import_module(args.deployment_script.resolve())
    registry = BundleRegistry(args.model_dir.resolve())
    rna, _, metadata = module.load_scp()
    rows_out: list[dict] = []

    for route in module.ROUTES:
        routed = module.routed_meta(metadata, route)
        external_rows = routed["external_row"].to_numpy(dtype=np.int64)
        positions = np.unique(
            np.linspace(0, len(external_rows) - 1, min(args.sample_cells, len(external_rows)), dtype=np.int64)
        )
        sample_rows = external_rows[positions]
        for modality in module.MODALITIES:
            bundle = registry.bundle(route["held_out_cell_type"], modality)
            observed_latent, observed_distance = predict_latent(
                rna[sample_rows, :], bundle, n_jobs=args.n_jobs
            )
            original_dir = module.OUT / modality / route["route_name"]
            expected_latent = np.load(
                original_dir / "virtual_modality_latent.float32.npy",
                mmap_mode="r",
                allow_pickle=False,
            )[positions]
            expected_distance = np.load(
                original_dir / "rna_latent_knn25_mean_distance.float32.npy",
                mmap_mode="r",
                allow_pickle=False,
            )[positions]
            latent_error = np.abs(observed_latent - expected_latent)
            distance_error = np.abs(observed_distance - expected_distance)
            max_latent_error = float(np.max(latent_error))
            max_distance_error = float(np.max(distance_error))
            passed = max(max_latent_error, max_distance_error) <= args.tolerance
            rows_out.append(
                {
                    "cell_type": route["held_out_cell_type"],
                    "route_name": route["route_name"],
                    "modality": modality,
                    "n_sample_cells": int(len(positions)),
                    "max_abs_latent_error": max_latent_error,
                    "mean_abs_latent_error": float(np.mean(latent_error)),
                    "max_abs_knn_distance_error": max_distance_error,
                    "tolerance": args.tolerance,
                    "status": "PASS" if passed else "FAIL",
                }
            )
            print(f"{route['route_name']:12s} {modality:9s} {'PASS' if passed else 'FAIL'}")

    results = pd.DataFrame(rows_out)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, sep="\t", index=False)
    summary = {
        "status": "PASS" if (results["status"] == "PASS").all() and len(results) == 30 else "FAIL",
        "expected_comparisons": 30,
        "observed_comparisons": int(len(results)),
        "passed_comparisons": int((results["status"] == "PASS").sum()),
        "validation_table": str(args.output),
    }
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
