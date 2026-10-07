# DearDTI reuse in the fresh native Kd benchmark

Source: https://github.com/kangk1204/DearDTI

Frozen commit: `cc741b8adfa445e45dedfcef884567cb125a819a` (local checkout verified).

The selected GINE and amino-acid CNN encoders derive from
`experiments/plm_model.py` (original SHA256
`a989b656447be4f06c50f80779c60b38378153370a0d8dca490dd12ab54f763c`).
Cross-attention and bilinear fusion are copied from
`dti_sota/models/fusion.py` (SHA256
`68743fc5c550828c8a74feccd19487bf2e4db56f29e246fc73bddc73fc9709ee`).

The full MIT notice is preserved in
`src/nutriomics_dti/vendor/DEARDTI_MIT_LICENSE.txt`.

Changes: GINE BatchNorm excludes padded atoms and handles a singleton valid atom;
CNN intermediate states exclude padded residues; adaptive pooling uses each real
sequence's boundaries. New RDKit atom/bond features, exact native Kd labels and one
pKd regression head replace the old data interfaces and multi-task heads. Full
protein sequences identify targets; the sequence encoder uses the first 1,024
residues. Molecules above 128 atoms are explicitly excluded before all benchmark
splits, equally for controls and GINE, and counted. No molecular graph is silently
truncated. The fresh amino-acid encoder uses no pretrained language model.

No old weights, KIBA transformations, original test-selected recipes, affinity
scores, retrieval caches or uncertainty calibration are reused. Validation-only
parameter/epoch selection precedes one final test attempt for each declared model.
Current-release internal holdouts do not establish independence from historical
models whose complete source ledger remains unavailable.
