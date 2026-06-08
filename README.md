# Catalan_Mg

This repository now has a reproducible generation pipeline for the exploratory Catalan / Mg / ZDED visuals that were previously stored mostly as flat artifacts.

The PDF note at `catalan_ruelle_mayer_research_note_final_with_mg_figures.pdf` is still the narrative source. The code under `src/catalan_mg/` now rebuilds the figure families, data tables, and SVG diagnostics into an organized `results/` tree.

## Layout

```text
.
|- catalan_ruelle_mayer_research_note_final_with_mg_figures.pdf
|- src/catalan_mg/                 # canonical generators
|- generate_results.py             # main entrypoint
|- MG_lattice.py                   # legacy compatibility entrypoint
|- zded_mg_lattice_calibration.py  # legacy compatibility entrypoint
|- results/
|  |- figures/
|  |- data/
|  |- html/
|  `- generated_manifest.csv
`- requirements.txt
```

The historical root-level PNG/SVG/CSV files are still present for compatibility. The canonical output location is now `results/`.

## Rebuild

Install dependencies if needed:

```bash
python3 -m pip install -r requirements.txt
```

Generate the organized results tree:

```bash
python3 generate_results.py --output-root results --layout canonical
```

Recreate the historical flat-layout outputs for the original scripts:

```bash
python3 MG_lattice.py
python3 zded_mg_lattice_calibration.py
```

## Generated outputs

- `results/figures/` contains the MG lattice plots, Catalan radical board, prime diagnostics, pressure plots, and the exploratory ZDED toy calibration figure.
- `results/data/` contains the MG lattice CSV exports and the deterministic ZDED irrational-layer hypothesis scan.
- `results/html/` contains the Lyapunov-vs-Catalan SVG variants.
- `results/generated_manifest.csv` lists every artifact produced by the canonical build.

## Scope guardrails

The PDF explicitly treats these figures as exploratory debugging surfaces, not proof objects. This repo now follows that same boundary:

- Safe: reproducible Mg-channel visualizations, pressure proxy scans, Catalan radical convergence boards, Ti2/tanh prime diagnostics, and toy calibration experiments.
- Unsafe: claiming a direct Selberg-zeta identity for Catalan's constant, treating Liouville numbers as singularities of a signed determinant, or presenting the ZDED scan as an exact collapse/proof.

## Notes

- The pressure plots implement a finite-cutoff proxy for the branch-summability threshold `u = 2t - 1` discussed in the PDF.
- The `Ti2` plots use the dilogarithm relation `Ti2(x) = Im Li_2(ix)`.
- The ZDED plot is deliberately labeled as a toy calibration because the current model remains exploratory and does not tightly collapse the rational layers.
