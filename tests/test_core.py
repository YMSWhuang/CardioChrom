import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from cardiochrom.bundle import FrozenBundle
from cardiochrom.core import decode_features, log_normalize, predict_latent


def make_bundle(tmp_path, modality="ATAC"):
    rng = np.random.default_rng(7)
    path = tmp_path / modality
    path.mkdir()
    arrays = {
        "rna_components": rng.normal(size=(50, 36_100)).astype(np.float32),
        "train_rna_latent": rng.normal(size=(30, 50)).astype(np.float32),
        "train_target_latent": rng.normal(size=(30, 50)).astype(np.float32),
        "target_components": rng.normal(size=(50, 7)).astype(np.float32),
    }
    for name, value in arrays.items():
        np.save(path / f"{name}.npy", value, allow_pickle=False)
    (path / "metadata.json").write_text(
        json.dumps({"cell_type": "Endothelial", "modality": modality, "knn_k": 25}),
        encoding="utf-8",
    )
    return FrozenBundle.load(path), arrays


class CoreTests(unittest.TestCase):
    def test_sparse_log_normalization(self):
        counts = sp.csr_matrix(np.asarray([[1, 1, 0], [0, 0, 0]], dtype=np.float32))
        observed = log_normalize(counts).toarray()
        expected = np.asarray([[np.log1p(5000), np.log1p(5000), 0], [0, 0, 0]], dtype=np.float32)
        np.testing.assert_allclose(observed, expected, rtol=1e-6)

    def test_predict_latent_shapes_and_finiteness(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle, _ = make_bundle(Path(directory))
            counts = sp.csr_matrix((np.ones(4), ([0, 0, 1, 1], [0, 10, 3, 20])), shape=(2, 36_100))
            latent, distance = predict_latent(counts, bundle)
            self.assertEqual(latent.shape, (2, 50))
            self.assertEqual(distance.shape, (2,))
            self.assertTrue(np.isfinite(latent).all())
            self.assertTrue(np.isfinite(distance).all())

    def test_atac_decode_is_linear(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle, arrays = make_bundle(Path(directory), "ATAC")
            latent = np.ones((2, 50), dtype=np.float32)
            observed = decode_features(latent, bundle, [1, 4])
            expected = latent @ arrays["target_components"][:, [1, 4]]
            np.testing.assert_allclose(observed, expected, rtol=1e-6, atol=1e-6)

    def test_histone_subset_uses_full_denominator(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle, arrays = make_bundle(Path(directory), "H3K27ac")
            latent = np.ones((1, 50), dtype=np.float32)
            raw = latent @ arrays["target_components"]
            raw = np.expm1(np.clip(raw, 0, 20))
            expected_full = np.log1p(raw * (10_000 / raw.sum(axis=1, keepdims=True)))
            observed = decode_features(latent, bundle, [0, 3, 6])
            np.testing.assert_allclose(observed, expected_full[:, [0, 3, 6]], rtol=1e-5, atol=1e-5)


if __name__ == "__main__":
    unittest.main()
