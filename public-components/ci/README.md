# Bounded public-component CI

The `Public components` GitHub Actions workflow runs on a standard Ubuntu CPU
runner with Python 3.12 and Node 22. It runs for changes to this distribution or
the workflow, and can be started manually. The exact run and commit are the
source of truth for remote success; `VALIDATION.json` records local checks.

The workflow verifies the exact Git-tracked release inventory and file hashes,
then runs the GET/HEAD API tests, integration contracts, food normalization,
food dictionary and five named evidence-parser test files. The default
integration dispatcher is exercised through a tiny synthetic SQLite query in
both public and original directory layouts, including conflicting import
environments. Results, actual first-party module paths/hashes and input
immutability are checked.

The methylation export and intervention validation regressions compile their
exact native function ASTs and use native file I/O helpers. They exercise small
synthetic output tables without importing model/scientific dependencies. These
checks do not establish full package importability, model performance, real
study validity or clinical utility. Native BH logic and analysis protocols are
unchanged; explicitly non-estimable rows remain supported.

Export regressions also preserve seeded obsolete files and removed-input
leftovers while excluding them from the current-call receipt. Atlas regressions
distinguish missing or filename-only inventories from invalid and hash-verified
receipts. CLI tests use the native entrypoint with the exact isolated validator
bound into its import slot; they verify exit codes 2, 1 and 0 without importing
SciPy. This does not establish the full scientific package runtime.

`requirements.txt` pins the nine direct dependencies used locally. Transitive
dependencies are not locked. Installation accepts wheels only. No scientific
package extras, model frameworks, data downloads, training, caches or uploaded
artifacts are part of this workflow. The public API needs no npm installation.
Action revisions are pinned, checkout does not persist credentials, and the
job has read-only contents permission and a ten-minute limit.

The four R workers and the remaining scientific test suites are outside this
CI gate. Scientific failures, historical result tables, submission readiness
and deployment status retain their separate evidence requirements.
