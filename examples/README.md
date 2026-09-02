# Minimal synthetic example

This example checks that the installed software, extracted model bundle, input
format, and command-line interface work together. It creates three synthetic
cells (`vCM`, `Fibroblast`, and `Endothelial`) with sparse integer counts across
the exact frozen gene list.

The generated values are random and are **not a biological dataset**. Do not use
the predictions for scientific interpretation.

## 1. Install CardioChrom

From the repository root, follow the full environment instructions in the main
README. The standard installation commands are:

```bash
python -m pip install .
cardiochrom --help
```

On shared systems, unrelated user-site packages can leak into an otherwise new
environment. If `python -m pip check` reports packages that CardioChrom did not
install, isolate the environment before installation:

```bash
unset PYTHONPATH
export PYTHONNOUSERSITE=1
```

## 2. Extract the model archive

Download the versioned model archive and verify its checksum before extraction.
After extraction, the directory passed to `--model-dir` should contain
`common/genes.txt`, `bundle_manifest.json`, and the fold-specific bundles.

## 3. Create a tiny input file

```bash
python examples/create_toy_input.py \
  --model-dir /path/to/CardioChrom_model_bundle_v1 \
  --output example_input.h5ad
```

The generated `example_input.h5ad` contains raw non-negative integer counts and
a `cell_type` column in `obs`.

## 4. Run frozen inference

```bash
cardiochrom predict \
  --input example_input.h5ad \
  --model-dir /path/to/CardioChrom_model_bundle_v1 \
  --cell-type-key cell_type \
  --modalities ATAC,H3K27ac,H3K27me3 \
  --output example_results
```

Because the example uses canonical CardioChrom cell-type names, no routing map is
needed. Real datasets with different labels must provide an explicit two-column
routing map with `input_label` and `canonical_cell_type`.

## 5. Check the outputs

A successful run creates:

- `example_results/CardioChrom_latent_predictions.npz`: 50-dimensional
  predictions and KNN25 mean distances for each requested modality;
- `example_results/CardioChrom_cells.tsv.gz`: input labels, canonical routes,
  and routing status for each cell;
- `example_results/CardioChrom_run_manifest.json`: input coverage, route
  counts, modalities, normalization, and translator metadata.

The validated smoke test produces 3 routed cells, frozen-gene coverage of 1.0,
no missing frozen genes, latent arrays with shape `(3, 50)`, distance arrays
with shape `(3,)`, and finite values in every output array.

For real data, supply an `.h5ad` containing raw RNA counts, unique gene symbols,
and a cell-type column. See the root README for the full input contract.
