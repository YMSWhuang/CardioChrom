# Public-release checklist

- [x] Export the frozen portable model bundle.
- [x] Validate the portable bundle against the frozen reference outputs (30/30 comparisons passed).
- [x] Review `checksums.sha256`, `bundle_manifest.json`, and the validation table.
- [x] Remove internal filesystem paths and infrastructure-specific scripts from the GitHub repository.
- [x] Add a minimal synthetic example-input generator and quick-start instructions.
- [x] Run the documented end-to-end quick start in a clean Python 3.11 environment (3/3 cells routed, gene coverage 1.0, all outputs finite).
- [x] Make the source repository publicly accessible.
- [x] Add the MIT open-source license and citation metadata.
- [x] Reserve Zenodo DOI `10.5281/zenodo.22239579` for the portable model bundle and add it to repository metadata.
- [ ] Publish the model archive under the reserved Zenodo DOI.
- [ ] Add the final paper citation and Data Availability Statement.
- [ ] Tag the submission version as `v1.0.0-paper`.
