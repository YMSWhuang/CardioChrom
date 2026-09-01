# Portable model bundle

The model archive is separate from the GitHub code repository because the exact
frozen assets are approximately 3 GB before archive compression.

```text
CardioChrom_model_bundle_v1/
  bundle_manifest.json
  checksums.sha256
  common/
    genes.txt
    atac_features.tsv.gz
    histone_bins.tsv.gz
  folds/
    00_vCM/
      ATAC/
        metadata.json
        rna_components.npy
        train_rna_latent.npy
        train_target_latent.npy
        target_components.npy
      H3K27ac/
        ...
      H3K27me3/
        ...
    ...
    11_Neuronal/
```

Each modality has the same inference contract:

1. align raw RNA counts to `common/genes.txt`;
2. normalize each cell to 10,000 and apply `log1p`;
3. project with `rna_components.npy`;
4. apply uniform Euclidean KNN25 using `train_rna_latent.npy` as predictors and
   `train_target_latent.npy` as targets;
5. retain the predicted 50-dimensional target latent and the mean KNN25 distance;
6. optionally decode with `target_components.npy`.

H3K27me3 is made portable by exporting the reconstructed train-only RNA and target
latents. Public inference must not depend on the original development RNA matrix,
H3K27me3 matrix, modeling manifest, donor labels, or split metadata.

