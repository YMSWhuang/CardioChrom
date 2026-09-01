# Public-release checklist

- [ ] Run `scripts/export_frozen_bundles_on_puma.py`.
- [ ] Run `scripts/validate_against_scp3342_on_puma.py` and require all 30 route-by-modality comparisons to pass.
- [ ] Review `checksums.sha256`, `bundle_manifest.json`, and the validation TSV.
- [ ] Add a very small redistributable example RNA dataset.
- [ ] Run the documented end-to-end quick start in a clean environment.
- [ ] Confirm University of Arizona IP/public-disclosure timing.
- [ ] Add the approved open-source license.
- [ ] Create a private GitHub repository and remove internal filesystem paths from release-facing files.
- [ ] Deposit the model archive and frozen paper version in Zenodo.
- [ ] Insert the model DOI and paper-code DOI in the README and Data Availability Statement.
- [ ] Tag the submission version as `v1.0.0-paper`.

