# Plan

This plan is based on the PDF note in the repo `docs`.

The working assumptions are:

- zero proof inflation
- explicit labeling of speculative numerics
- scoped exploratory generators instead of narrative-only artifacts
- reproducibility first, claims second

## Current objective

Turn the repo into a reproducible exploratory lab for the PDF's three main tracks:

1. the signed / twisted Catalan determinant track
2. the separate positive pressure / cusp track
3. the MG-lattice visual and numerical scaffolding

## Completed

- Added a canonical source tree under `src/catalan_mg/`.
- Added `scripts/generate_results.py` to rebuild the organized `results/` tree.
- Fixed `REPO_ROOT` path bug in `generate_results.py` (was pointing at `scripts/` instead of repo root).
- Removed stale `results/zded_irrational_hypotheses.csv` duplicate (canonical copy lives in `results/data/`).
- Kept `scripts/MG_lattice.py` and `scripts/zded_mg_lattice_calibration.py` as compatibility entry-points.
- Rebuilt missing generator coverage for the orphan PNG/SVG artifacts.
- Added `README.md` and `RESULTS.md`.
- All 24 artifacts (14 figures, 3 CSVs, 6 SVGs, manifest) now regenerate cleanly from one command:
  ```
  python scripts/generate_results.py
  ```

## Next workstreams

### 1. Determinant reproducibility

- Implement the beta-matrix determinant experiments from the PDF appendix with stable precision controls.
- Sweep over matrix size `N`, precision, and parameter `t`; write outputs to `results/data/`.
- Add one comparison plot showing numerical stability across `N` rather than a single snapshot.
- Use `mpmath` for high-precision arithmetic; record the working precision alongside each result.

### 2. Trace and periodic-word checks

- Enumerate periodic odd-continued-fraction words for short lengths `m = 2, 3`.
- Verify sign character bookkeeping `chi(w) = (-1)^(sum k_i)` numerically against the PDF's claim.
- Export trace tables to `results/data/` so they can be inspected independently of any plot.

### 3. Pressure track

- Extend the current finite-cutoff proxy from a single sweep into a small grid over cutoff values.
- Detect where the pressure picture is stable versus where it is a truncation artifact.
- Keep this track in its own section; do not mix with signed-determinant claims.
- Report the cutoff sensitivity explicitly in the figure titles.

### 4. ZDED track

- Treat the current ZDED output as a toy calibration benchmark only.
- Improve the parameter search (larger grid or simple gradient descent) and write residuals for both rational and irrational layers.
- Do not promote this track beyond "exploratory model comparison" until rational layers collapse tightly.
- Add a second CSV column for the rational-layer calibration residuals so the calibration quality is machine-checkable.

### 5. Test coverage

- Add a lightweight smoke test (e.g. `pytest tests/`) that runs `generate_all_artifacts` on a small `--section mg` call and asserts the manifest is non-empty and all listed files exist on disk.
- This catches path bugs and import regressions without re-running the full 2-minute generation.

### 6. Proof hygiene

- Maintain a clear boundary between "safe statement", "conjectural statement", and "visual heuristic".
- When new claims are added, cite the exact PDF section or numerical artifact that supports them.
- Keep the quarantine list from the PDF inside `docs/`.

## Success criteria

- A fresh clone can rebuild all figures and CSVs with one documented command.
- MG-lattice outputs and Catalan radical outputs are stable across reruns.
- Determinant, trace, pressure, and ZDED experiments are stored as separate reproducible tracks under `results/`.
- Speculative claims stay labeled as speculative.

## Do not

- Claim a direct Selberg zeta identity for Catalan's constant from the current numerics.
- Treat Liouville numbers as poles or singularities of the signed Catalan determinant.
- Treat fitted exponents or finite-cutoff visual agreement as proof.
- Describe the current ZDED run as an exact calibration.
