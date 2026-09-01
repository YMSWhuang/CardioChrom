from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import KNeighborsRegressor

from .bundle import FrozenBundle
from .constants import KNN_K, NORMALIZATION_SCALE


def log_normalize(
    counts: sp.spmatrix | np.ndarray,
    library_size: np.ndarray | None = None,
    scale: float = NORMALIZATION_SCALE,
) -> sp.csr_matrix | np.ndarray:
    if sp.issparse(counts):
        matrix = counts.tocsr().astype(np.float32, copy=True)
        if matrix.data.size and (not np.isfinite(matrix.data).all() or np.min(matrix.data) < 0):
            raise ValueError("RNA counts must be finite and non-negative")
        lib = np.asarray(matrix.sum(axis=1)).ravel() if library_size is None else np.asarray(library_size).ravel()
        if lib.size != matrix.shape[0]:
            raise ValueError("Library-size length does not match the number of cells")
        factors = np.zeros_like(lib, dtype=np.float64)
        positive = lib > 0
        factors[positive] = scale / lib[positive]
        matrix = matrix.multiply(factors[:, None]).tocsr()
        matrix.data = np.log1p(matrix.data).astype(np.float32, copy=False)
        return matrix

    matrix = np.asarray(counts, dtype=np.float32)
    if matrix.ndim != 2 or not np.isfinite(matrix).all() or np.min(matrix) < 0:
        raise ValueError("RNA counts must be a finite, non-negative two-dimensional matrix")
    lib = matrix.sum(axis=1, dtype=np.float64) if library_size is None else np.asarray(library_size).ravel()
    if lib.size != matrix.shape[0]:
        raise ValueError("Library-size length does not match the number of cells")
    factors = np.zeros_like(lib, dtype=np.float64)
    positive = lib > 0
    factors[positive] = scale / lib[positive]
    return np.log1p(matrix * factors[:, None]).astype(np.float32)


def predict_latent(
    counts: sp.spmatrix | np.ndarray,
    bundle: FrozenBundle,
    n_jobs: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    normalized = log_normalize(counts)
    rna_components = np.asarray(bundle.array("rna_components"), dtype=np.float32)
    query_latent = np.asarray(normalized @ rna_components.T, dtype=np.float32)
    train_rna = np.asarray(bundle.array("train_rna_latent"), dtype=np.float32)
    train_target = np.asarray(bundle.array("train_target_latent"), dtype=np.float32)
    knn = KNeighborsRegressor(
        n_neighbors=KNN_K,
        weights="uniform",
        metric="euclidean",
        algorithm="brute",
        n_jobs=max(1, int(n_jobs)),
    )
    knn.fit(train_rna, train_target)
    prediction = knn.predict(query_latent).astype(np.float32, copy=False)
    distances, _ = knn.kneighbors(query_latent, n_neighbors=KNN_K, return_distance=True)
    mean_distance = distances.mean(axis=1).astype(np.float32)
    return prediction, mean_distance


def decode_features(
    latent: np.ndarray,
    bundle: FrozenBundle,
    feature_indices: np.ndarray | list[int] | None = None,
    chunk_size: int = 256,
) -> np.ndarray:
    latent = np.asarray(latent, dtype=np.float32)
    if latent.ndim != 2 or latent.shape[1] != 50:
        raise ValueError("Latent matrix must have shape (cells, 50)")
    components = bundle.array("target_components")
    indices = None if feature_indices is None else np.asarray(feature_indices, dtype=np.int64)
    if indices is not None:
        if indices.ndim != 1 or np.any(indices < 0) or np.any(indices >= components.shape[1]):
            raise ValueError("Feature indices are out of range")

    if bundle.modality == "ATAC":
        decoder = np.asarray(components if indices is None else components[:, indices], dtype=np.float32)
        return np.asarray(latent @ decoder, dtype=np.float32)

    output_width = components.shape[1] if indices is None else len(indices)
    output = np.empty((latent.shape[0], output_width), dtype=np.float32)
    decoder = np.asarray(components, dtype=np.float32)
    for start in range(0, latent.shape[0], chunk_size):
        end = min(start + chunk_size, latent.shape[0])
        reconstructed = np.asarray(latent[start:end] @ decoder, dtype=np.float32)
        np.clip(reconstructed, 0.0, 20.0, out=reconstructed)
        np.expm1(reconstructed, out=reconstructed)
        denominator = reconstructed.sum(axis=1, keepdims=True, dtype=np.float64)
        denominator = np.maximum(denominator, 1e-12)
        reconstructed *= (NORMALIZATION_SCALE / denominator).astype(np.float32)
        np.log1p(reconstructed, out=reconstructed)
        output[start:end] = reconstructed if indices is None else reconstructed[:, indices]
    return output

