from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .constants import CANONICAL_CELL_TYPES, KNN_K, LATENT_DIM, MODALITIES, N_GENES


@dataclass(frozen=True)
class FrozenBundle:
    path: Path
    cell_type: str
    modality: str
    metadata: dict

    @classmethod
    def load(cls, path: str | Path) -> "FrozenBundle":
        bundle_path = Path(path)
        metadata_path = bundle_path / "metadata.json"
        if not metadata_path.is_file():
            raise FileNotFoundError(metadata_path)
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        bundle = cls(
            path=bundle_path,
            cell_type=str(metadata["cell_type"]),
            modality=str(metadata["modality"]),
            metadata=metadata,
        )
        bundle.validate()
        return bundle

    def array(self, name: str, mmap_mode: str | None = "r") -> np.ndarray:
        path = self.path / f"{name}.npy"
        if not path.is_file():
            raise FileNotFoundError(path)
        return np.load(path, mmap_mode=mmap_mode, allow_pickle=False)

    def validate(self) -> None:
        if self.cell_type not in CANONICAL_CELL_TYPES:
            raise ValueError(f"Unsupported canonical cell type: {self.cell_type}")
        if self.modality not in MODALITIES:
            raise ValueError(f"Unsupported modality: {self.modality}")
        expected = {
            "rna_components": (LATENT_DIM, N_GENES),
            "train_rna_latent": (None, LATENT_DIM),
            "train_target_latent": (None, LATENT_DIM),
            "target_components": (LATENT_DIM, None),
        }
        observed_rows = None
        for name, shape in expected.items():
            array = self.array(name)
            if array.ndim != 2:
                raise ValueError(f"{self.path}/{name}.npy must be two-dimensional")
            if shape[0] is not None and array.shape[0] != shape[0]:
                raise ValueError(f"{name} shape {array.shape} does not match {shape}")
            if shape[1] is not None and array.shape[1] != shape[1]:
                raise ValueError(f"{name} shape {array.shape} does not match {shape}")
            if name in {"train_rna_latent", "train_target_latent"}:
                if observed_rows is None:
                    observed_rows = array.shape[0]
                elif array.shape[0] != observed_rows:
                    raise ValueError("Train RNA and target latent row counts differ")
        if int(self.metadata.get("knn_k", -1)) != KNN_K:
            raise ValueError("Bundle does not declare the frozen KNN25 translator")


class BundleRegistry:
    def __init__(self, model_dir: str | Path):
        self.model_dir = Path(model_dir)
        manifest_path = self.model_dir / "bundle_manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(manifest_path)
        self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.genes = tuple(
            line.strip()
            for line in (self.model_dir / "common" / "genes.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
        if len(self.genes) != N_GENES or len(set(self.genes)) != N_GENES:
            raise ValueError("Frozen gene order must contain 36,100 unique genes")
        self._bundles: dict[tuple[str, str], FrozenBundle] = {}

    def bundle(self, cell_type: str, modality: str) -> FrozenBundle:
        key = (cell_type, modality)
        if key not in self._bundles:
            if cell_type not in CANONICAL_CELL_TYPES or modality not in MODALITIES:
                raise ValueError(f"Unsupported route: {cell_type}/{modality}")
            fold = CANONICAL_CELL_TYPES.index(cell_type)
            path = self.model_dir / "folds" / f"{fold:02d}_{cell_type}" / modality
            self._bundles[key] = FrozenBundle.load(path)
        return self._bundles[key]

