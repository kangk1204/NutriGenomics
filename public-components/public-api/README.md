# Public projection HTTP API

Entry point: `worker.mjs`. The handler imports only packaged public JSON snapshots. It accepts GET and HEAD and has no filesystem, queue, training, model-download or arbitrary artifact route.

| Route | Response |
|---|---|
| `/api/v1/health` | Actual release counts and source-manifest hash |
| `/api/v1/components` | Separate component origins and unverified model HTTP status |
| `/api/v1/provenance` | Source releases and public distribution scope |
| `/api/v1/compounds?limit=20&offset=0&q=` | Paginated complete structure identity and field provenance |
| `/api/v1/foods?limit=20&offset=0&q=` | Paginated native USDA food definitions |
| `/api/v1/foods/{encoded_id}` | One native food definition |
| `/api/v1/foods/{encoded_id}/reports?limit=20&offset=0` | Native quantities, units, basis and source records |
| `/api/v1/exports/compounds` | Exact approved compound snapshot |
| `/api/v1/exports/foods` | Derived pilot snapshot with upstream dataset hash and native measure definitions |

List responses have `schema_version`, `release_id`, `total`, `offset`, `limit`, and `items`. Errors have `schema_version` and `error.code/message`. Limit is a strict integer 1–100, offset 0–1,000,000. Repeated/unknown parameters, oversized search, nonintegers and invalid ranges return 400. Missing records/routes return 404. Mutations return 405 with `Allow: GET, HEAD`. Search produces candidates; it does not infer identity equivalence.

Run `npm test` with Node 22+; there are no npm dependencies. Tests use the packaged real cleared snapshots and verify success, error semantics, pagination, provenance, raw export and absence of worker controls. See `openapi.json` for the request/response contract.

No public API URL or new deployed commit is asserted by this package. To integrate with existing hosting, the exact Site owner must build this entrypoint into its configured Worker, preserve its project ID, push the exact source, save/deploy through the authorized route, and then verify unauthenticated HTTP requests against that deployed commit. Credential renewal is a separately pending action; this package creates no token, host trust or network configuration.
