# Distribution notice

The owner-authorized, independently authored additions in this directory are intentionally contributed under Apache-2.0; the complete license is in `LICENSE`. This does not claim that private upstream repositories previously had an Apache or MIT license. Preserve existing third-party and component-specific terms:

- `methylation/` retains its full MIT LICENSE, copyright 2026 NutriOmics contributors.
- DearDTI vendor encoders and fusion retain the full `compound-target/src/nutriomics_dti/vendor/DEARDTI_MIT_LICENSE.txt` and `NOTICE-DEARDTI.md`, copyright 2025 DearDTI contributors.
- Data retain their separate field-scoped USDA/NCBI rights manifests. ChEMBL attribution is preserved as a source-adapter notice; this release includes no ChEMBL data. Apache is not a data license for mixed external sources.

Changes from frozen upstream code are marked by distinct hashes in `SOURCE_MANIFEST.json`: strict artifact receipt and typed acceptance checks, graceful optional image model initialization, queue process ownership and atomic report revision, neural receipt drift detection, relation-generation/abstention separation, metric input guards, strict SQL/transaction behavior, exact dictionary mapping counts, and a corrected R contrast export label. Existing scientific results and original private branches are not rewritten.

Public packaging changes require an explicit R input root instead of a personal host path; remove the withheld Kang license package-data reference; correct the integration test's sibling package path and declare its local test imports; and omit historical script-dependent or private rebuild tests whose inputs are not part of this release. No scientific fitting or protocol change was performed.

KangDTI-derived helpers, control models, dependent CLI/tests and settings are withheld because their existing notice authorizes private research reuse only. Participant-level data/predictions, clinical metadata, reports, source literature bodies, unresolved FooDB/CTD/BindingDB/ChEMBL data, private seeds and trained model assets are withheld. Data-access application approval, including RDA academic/education use, does not by itself authorize public redistribution.

The original nine files in the public destination are preserved. This release imports no private Git history. Its source provenance identifies the distinct original components and commits; repository ownership/permissions do not establish a live runtime owner or deployment.
