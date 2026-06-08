# Results

Canonical build command:

```bash
python3 generate_results.py --output-root results --layout canonical
```

The current generated artifact set is indexed in `results/generated_manifest.csv`.

## Output summary

- MG lattice family:
  - `results/figures/mg_lattice_checkerboard.png`
  - `results/figures/mg_lattice_double_helix.png`
  - `results/figures/mg_lattice_prime_style_funnel.png`
  - `results/figures/catalan_apery_3d_heatmap.png`
  - `results/figures/catalan_Mg_3d.png`
  - `results/figures/Mg_catalan_zeta.png`
  - `results/data/mg_lattice_discrete.csv`
  - `results/data/mg_lattice_continuous.csv`

- Catalan radical / Ramanujan family:
  - `results/figures/catalan_radical_riemann.png`
  - `results/figures/ramanujan_catalan_spiral_projection.svg`
  - These correspond to the PDF's Catalan radical / Riemann-sphere convergence board theme.

- Pressure / cusp proxy family:
  - `results/figures/pressure_phase_diagram.png`
  - `results/figures/pressure_slice_t1.png`
  - These implement the positive pressure proxy around the branch threshold `u = 2t - 1`.

- Prime diagnostics:
  - `results/figures/Figure_1.png`
  - `results/figures/tanh_diophantine.png`
  - `results/figures/ti2_primes.png`

- ZDED calibration:
  - `results/figures/zded_mg_lattice_calibration.png`
  - `results/data/zded_irrational_hypotheses.csv`
  - Deterministic default run: the toy calibration residual is about `6.81e-01`, which is still exploratory rather than collapsed.

- HTML / SVG diagnostics:
  - `results/html/flight_map_lyapunov_vs_catalan.svg`
  - `results/html/flight_map_lyapunov_vs_catalan_bijective.svg`
  - `results/html/flight_map_lyapunov_vs_catalan_bijective_transposed.svg`
  - `results/html/flight_map_lyapunov_vs_catalan_injective.svg`
  - `results/html/flight_map_lyapunov_vs_catalan_injective_transposed.svg`
  - `results/html/flight_map_lyapunov_vs_catalan_surjective.svg`

## Interpretation

- The MG lattice figures are the clearest reproducible outputs tied directly to the PDF. They expose the parity switch between the rational and deep channels.
- The Catalan radical board is a visualization layer for convergence behavior, not a proof of arithmetic structure.
- The pressure figures are finite-cutoff scans meant to support the PDF's separate positive-pressure track for cusp behavior.
- The ZDED plot is intentionally scoped as a toy model. It is useful as a reproducible experiment, but not as evidence of a solved matching problem.
