# Plan

This plan is based on the PDF note in the repo  `docs`.

 The working assumption here is:

- zero proof inflation
- explicit labeling of speculative numerics
- scoped exploratory generators instead of narrative-only artifacts
- reproducibility first, claims second

## Current objective

Turn the repo into a reproducible exploratory lab for the PDF's three main tracks:

1. the signed / twisted Catalan determinant track
2. the separate positive pressure / cusp track
3. the MG-lattice visual and numerical scaffolding

## Completed by now

- Added a canonical source tree under `src/catalan_mg/`.
- Added `generate_results.py` to rebuild the organized `results/` tree.
- Kept `MG_lattice.py` and `zded_mg_lattice_calibration.py` as compatibility entrypoints.
- Rebuilt the missing generator coverage for the orphan PNG/SVG artifacts.
- Added `README.md` and `RESULTS.md`.

## Next workstreams

### 1. Determinant reproducibility

- Implemented> the beta-matrix determinant experiments from the PDF appendix with stable precision controls.
- Determinant sweeps by matrix size `N`, precision, and parameter `t` into `results/data/`.
- We Added one comparison plot showing numerical stability instead of a single snapshot.


### 2. Trace and periodic-word checks

- Enumerate periodic odd-continued-fraction words for short lengths `m = 2, 3`.
- Verify the sign character bookkeeping `chi(w) = (-1)^(sum k_i)` numerically.
- Export trace tables that can be inspected independently of the plots.

### 3. Pressure track

- Extend the current finite-cutoff proxy into a small sweep over cutoffs, not a single proxy line.
- Detect where the pressure picture is stable versus where it is an artifact of truncation.
- Keep this track separate from the signed determinant claims, matching the PDF's split architecture.

### 4. ZDED track

- Treat the present ZDED output as a toy calibration benchmark only.
- Improve the parameter search and write out residuals on both rational and irrational layers.
- Do not promote this track beyond "exploratory model comparison" until rational layers collapse tightly.

### 6. Proof hygiene

- We will Keep a clear boundary between "safe statement", "conjectural statement", and "visual heuristic".
- When new claims are added, we will cite the exact source section or numerical artifact that supports them.
- We will Maintain the quarantine list from the PDF inside the docs.

## Success criteria

- A fresh clone can rebuild the figures and CSVs with one documented command.
- MG-lattice outputs and Catalan radical outputs remain stable across reruns.
- Determinant, collocation, pressure, and ZDED experiments are stored as separate reproducible tracks.
- Speculative claims stay labeled as speculative.

## We Do-not>

- Claim A direct Selberg zeta identity for Catalan's constant from the current numerics.
- Treat Liouville numbers as poles or singularities of the signed Catalan determinant.
- Treat fitted exponents or finite-cutoff visual agreement as proof.
- Describe the current ZDED run as an exact calibration.
