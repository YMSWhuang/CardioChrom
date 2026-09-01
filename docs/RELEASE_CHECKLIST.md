# Public-release checklist

- [x] Export the frozen portable model bundle.
- [x] Validate the portable bundle against the frozen reference outputs (30/30 comparisons passed).
- [x] Review `checksums.sha256`, `bundle_manifest.json`, and the validation table.
- [x] Remove internal filesystem paths and infrastructure-specific scripts from the GitHub repository.
- [ ] Add a very small redistributable example RNA dataset.
- [ ] Run the documented end-to-end quick start in a clean environment.
- [ ] Confirm the final public-release timing.
- [ ] Add the approved open-source license.
- [ ] Deposit the model archive and frozen paper version in Zenodo.
- [ ] Insert the model DOI and paper-code DOI in the README and Data Availability Statement.
- [ ] Tag the submission version as `v1.0.0-paper`.
