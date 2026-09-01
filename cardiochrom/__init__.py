"""CardioChrom frozen virtual-epigenome inference."""

from .bundle import BundleRegistry, FrozenBundle
from .core import decode_features, log_normalize, predict_latent

__all__ = [
    "BundleRegistry",
    "FrozenBundle",
    "decode_features",
    "log_normalize",
    "predict_latent",
]

__version__ = "0.1.0rc2"
