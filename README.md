# CardioChrom

CardioChrom is a frozen cardiovascular virtual-epigenome inference framework.
It maps single-cell or single-nucleus RNA count profiles to three reconstructed
regulatory representations:

- CardioChrom-ATAC
- CardioChrom-H3K27ac
- CardioChrom-H3K27me3

The framework uses a fixed 36,100-gene RNA interface, fold-specific 50-dimensional
RNA and target representations, and a frozen KNN25 latent-state translator. It does
not retrain or tune on the deployment cohort.

## Release-candidate status

This repository is a packaging candidate, not yet a versioned public release. The
Python inference interface and a minimal synthetic smoke test are included, while
the frozen model bundle is distributed separately. Before the paper release:

1. publish the model archive in a versioned repository such as Zenodo;
2. run the documented end-to-end quick start in a clean environment;
3. confirm the final software license and release timing.

## Input contract

The CLI expects an `.h5ad` file containing raw non-negative RNA counts.

- rows: cells or nuclei;
- columns: gene symbols;
- `obs`: a cell-type column supplied with `--cell-type-key`;
- `X`, or a counts layer supplied with `--layer`;
- exact canonical cell types by default, or an explicit two-column routing map.

The strict default requires all 36,100 frozen genes. A lower coverage threshold may
be requested explicitly, in which case absent genes are filled with zero and the
coverage is recorded.

## Installation

```bash
python -m pip install -e .
```

## Quick start

```bash
cardiochrom predict \
  --input example.h5ad \
  --model-dir /path/to/CardioChrom_model_bundle_v1 \
  --cell-type-key cell_type \
  --routing-map examples/routing_map_scp3342.tsv \
  --modalities ATAC,H3K27ac,H3K27me3 \
  --output results
```

For a complete small test that first creates a valid synthetic `.h5ad`, follow
the [minimal synthetic example](examples/README.md).

Primary outputs are the exact 50-dimensional virtual-modality latent states and
RNA-latent KNN25 mean distances. The Python API can additionally decode selected
ATAC or histone features. ATAC values are continuous reconstructed accessibility
scores, not probabilities. Histone decoding follows the frozen candidate-A rule:
latent reconstruction, clipping to 0-20, `expm1`, per-cell normalization to 10,000,
and `log1p`.

## Supported canonical cardiac cell types

`vCM`, `aCM`, `Adipocyte`, `Fibroblast`, `Endothelial`, `Endocardial`,
`Epicardial`, `Pericyte`, `Myeloid`, `SM`, `Lymphoid`, and `Neuronal`.

CardioChrom does not silently infer or change cell-type routing. Dataset-specific
labels must be mapped explicitly.

## Repository layout

```text
cardiochrom/                  Python inference package
examples/                    Synthetic smoke test and routing-map example
tests/                       Synthetic numerical tests
docs/MODEL_BUNDLE.md         Portable model format
docs/RELEASE_CHECKLIST.md    Steps required before public release
```

## Scientific scope

CardioChrom produces RNA-derived reconstructed regulatory states through a learned
cross-modal prior. Its outputs should not be described as information independent of
RNA or as causal regulatory measurements.
