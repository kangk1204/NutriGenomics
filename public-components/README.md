# Reviewed Nutri research components

This directory is an additive public distribution of reviewed code, synthetic tests, and selected redistribution-cleared data. It is independent of this repository's older pipeline at its root. The seven components remain distinct; their exact upstream repositories and main commits are in `COMPONENTS.json`. Standalone food companions have no verified Git origin and are identified as such.

| Directory | Included purpose | Public HTTP inference |
|---|---|---|
| integration | Research job contracts, provenance, image-candidate adapter and report-package utilities | Unverified |
| evidence | Molecular evidence extraction, receipt integrity and source adapters | Unverified |
| methylation | Leakage-controlled methylation research code | Unverified |
| intervention-atlas | Intervention analysis code and four R workers | Unverified |
| compound-target | Compound-target code, validated metric contracts and DearDTI encoders | Unverified |
| food-quality | Structural food-data validation and transactional normalization | Unverified |
| food-dictionary | Stable concept IDs, source-scoped mappings and dictionary evaluation | Unverified |

`public-api/` contains a separate GET/HEAD-only Worker over cleared public snapshots. Its local request-response tests pass. This distribution does not assert that it is deployed: `deployment_status` and model inference remain unverified. The operational integration FastAPI app manages private research jobs and is not the public service entrypoint.

The bounded public data consist of 53 government-computed PubChem structures, two existing USDA pilot foods, and 172 native food-nutrient reports. Those reports preserve amounts, units, measurement basis, provenance, and `molecule_presence_asserted=false`. The USDA selection is a pilot, not a representative population sample. Complete compound InChIKeys remain distinct. Two energy measures are separately typed from the 125 nutrient descriptors; no nutrient-to-molecule equivalence is fabricated.

The separate research quantity criterion is **4,630 = 2,563 + 2,067**, satisfied under the 2026 v1.3 §4 count definition. The strict 263 subset is additional. These counts do not mean that all underlying chemical payloads have redistribution permission. This public snapshot contains 53 chemical structures. No clinical accuracy or dietary efficacy claim follows from this release.

## Reproduce bounded checks

[Bounded CI](ci/README.md) documents the exact CPU-only test gate and its
limitations. The repository's `Public components` Actions run verifies tracked
release hashes and runs the selected synthetic checks for its exact commit.

Use an existing Node 22+ runtime, with no npm installation:

```sh
cd public-components/public-api
npm test
```

The integration tests require Python 3.12+, its declared test dependencies and the sibling evidence package. With dependencies already available:

```sh
cd public-components/integration
python -m pytest -q tests
cd ../food-quality
python -m pytest -q test_normalize_food_db.py
cd ../food-dictionary
python -m pytest -q test_dictionary.py
```

Install each scientific package from its own directory when needed; dependencies and optional model frameworks are declared in its `pyproject.toml`. Scientific training, model downloads, raw cohort acquisition and full private-data rebuilds are separate operations. They were not executed for this public release. Synthetic tests establish code contracts, not scientific performance. Script-dependent tests whose scripts were outside the reviewed allowlist and four private dictionary rebuild tests are omitted and recorded in the source manifest.

`SOURCE_MANIFEST.json` distinguishes frozen upstream hashes, corrected review-input hashes and actual distribution hashes. `PUBLIC_MANIFEST.json` hashes every published file except itself. `NOTICE.md` identifies licensing, reviewed changes and exclusions. Raw public snapshots and their rights manifests are available in `public-data/`; no private history, participant results, reports, weights, research DB or credentials are included.

The export receipt now inventories `SUMMARY.md` after writing it. Intervention
validation checks the recorded study-family BH values, estimable effects and
uncertainty, counts, identities and declared receipt inventories. Default
integration dispatch supports one unambiguous public or original source layout
and checks executed first-party module paths/hashes against its source snapshot.
These repairs do not change scientific estimators or frozen protocol IDs.

DearMeal has no verified repository/harness identity in this release. The image adapter is not treated as DearMeal or as another similarly named project. API deployment must be coordinated with the exact existing Site owner and verified with unauthenticated HTTP success/error/schema/export checks against the actual deployed commit.
