# Carbon monoxide dose-response analysis

Analysis code and data for the manuscript:

> **Sub-lethal metabolic suppression by carbon monoxide near occupational exposure
> limits: a hierarchical Bayesian dose-response analysis in human brain cells.**
> Megan Hird, Prachi Stafford, Alexander G. Fletcher, Mari Herigstad.

A single, self-contained Python script (`co_analysis.py`) reproduces every
quantitative figure and table in the paper from one data workbook, fitting two
hierarchical Bayesian dose-response models (metabolic activity from resazurin; cell
death from Annexin V/PI flow cytometry) and the metabolic–death decoupling.

## Repository contents

```
.
├── co_analysis.py              # the analysis (loaders, models, figures, table, CLI)
├── requirements.txt            # pip dependencies
├── environment.yml             # conda environment (recommended on Windows)
├── setup.cfg                   # flake8 style config (PEP 8, line length 99)
├── CITATION.cff                # how to cite this repository
├── LICENSE                     # MIT licence (code)
├── data/
│   ├── Paper_Data_0_3.xlsx         # consolidated workbook (main analysis input)
│   ├── Paper_Data_Resazurin.xlsx   # raw resazurin (for the provenance check)
│   ├── DATA_DICTIONARY.md          # description of every sheet and its columns
│   └── LICENSE                     # CC BY 4.0 licence (data)
└── outputs/
    └── table2_estimates.csv        # reference output, for checking reproduction
```

## What it produces

| Output | File | Needs model fit? |
|--------|------|------------------|
| Figure 2 — metabolic dose-response (per line) | `figure2_metabolic.png` | yes |
| Figure 3 — cell death / viability (flat) | `figure3_death.png` | yes |
| Figure 4 — metabolic–death decoupling | `figure4_decoupling.png` | yes |
| Figure 5 — protein normalisation | `figure5_protein.png` | yes |
| Table 2 — per-line posterior estimates | `table2_estimates.csv` | yes |
| Supplementary — convergence table (R-hat, ESS) | `supp_convergence.csv` | yes |
| Supplementary — trace plots (per model) | `supp_trace_*.png` | yes |
| Supplementary — posterior predictive checks | `supp_ppc_*.png` | yes |

Figure numbering matches the manuscript. **Figure 1** is the experimental-design /
analysis schematic, produced separately in a vector editor; it is not generated here.
The supplementary diagnostics are built with `--diagnostics` (see below).

Every figure (data and supplementary) is written three times from the same draw: a
300-dpi `.png` (raster preview / submission) and vector `.pdf` and `.svg` (for final
production and hand-editing, with fonts embedded). The table above lists the PNG
names; the PDF and SVG carry the same basenames.

## Install

With pip:

```bash
pip install -r requirements.txt
```

With conda (recommended on Windows, where it installs PyMC's compiler toolchain cleanly):

```bash
conda env create -f environment.yml
conda activate co-dose-response
```

## Run

From the repository root (the script reads `data/Paper_Data_0_3.xlsx` by default):

```bash
python co_analysis.py --all --data data/Paper_Data_0_3.xlsx --outdir outputs
```

That fits both models once (4 chains × 2000 draws — the production settings;
~5 min on a normal multi-core machine), caches the posteriors, and writes all
figures and the table to `outputs/`.

Rebuild individual pieces (reusing the cached fit, so these are near-instant):

```bash
python co_analysis.py --figure 2 --data data/Paper_Data_0_3.xlsx --outdir outputs
python co_analysis.py --figure 2 4 5 --data data/Paper_Data_0_3.xlsx --outdir outputs
python co_analysis.py --table --data data/Paper_Data_0_3.xlsx --outdir outputs
python co_analysis.py --refit --all --data data/Paper_Data_0_3.xlsx --outdir outputs   # force a fresh fit
python co_analysis.py --help
```

Options: `--data PATH`, `--outdir DIR`, `--figure N [N ...]` (2–5), `--table`,
`--all`, `--refit`, `--draws`, `--tune`, `--chains`, `--cores`, `--seed`,
`--check-resazurin`, `--raw-data PATH`.

### The fit cache

The first run that needs a fit samples both models and caches the posteriors as
`_cache_resazurin.nc` / `_cache_composition.nc` (plus `_cache_meta.json`) in the
output folder. Later runs that only rebuild a figure or the table reuse the cache
and finish in seconds. The cache is invalidated automatically — and the models
re-fitted — if you pass `--refit`, change the data file, or change the sampling
settings (`--draws/--tune/--chains/--seed`). These cache files are git-ignored.

### Resazurin provenance check (optional)

The analysis reads metabolic activity from the processed `(3) Resazurin` sheet of
the main workbook. To confirm that processing is reproducible from source,
`--check-resazurin` rebuilds the corrected values from the raw resazurin workbook
(corrected = raw − mean(media) − mean(blank)) and verifies them against the
workbook's own "Corrected R" column (no model fit required):

```bash
python co_analysis.py --check-resazurin --raw-data data/Paper_Data_Resazurin.xlsx
```

It reports the maximum discrepancy — under 1 unit for every well except one
documented contaminated plate (a constant ~182-unit residual reflecting a stored
blank differing slightly from the recorded one) — and lists any wells flagged
(red fill) and excluded.

### Supplementary diagnostics (optional)

`--diagnostics` writes the model-checking material for the supplementary information:
a per-parameter convergence table (R-hat, bulk-ESS, tail-ESS for both models),
trace plots of the population-level hyperparameters, and posterior predictive
checks (resazurin: observed vs replicated log-fluorescence density; composition:
observed vs predicted quadrant proportions). It uses the fit (or the cache):

```bash
python co_analysis.py --diagnostics --data data/Paper_Data_0_3.xlsx --outdir outputs
```

Outputs: `supp_convergence.csv`, `supp_trace_resazurin.png`,
`supp_trace_composition.png`, `supp_ppc_resazurin.png`, `supp_ppc_composition.png`.
Run it at the production settings (the defaults) so the diagnostics correspond to
the posterior of record.

## Reproducibility notes

* The production run used 4 chains × 2000 draws (after 2000 tuning iterations),
  `target_accept = 0.95`. Convergence: resazurin max R-hat 1.000, composition 1.010,
  no divergences. Credible intervals are 94% equal-tailed. A short smoke test
  (`--draws 200 --chains 2`) gives the same point estimates but higher R-hat and is
  not for publication.
* Exact package versions used for the published run are recorded in the manuscript
  Methods. `requirements.txt` / `environment.yml` pin the minimum major versions
  (PyMC ≥ 5). Report your own versions with
  `python -c "import pymc, arviz; print(pymc.__version__, arviz.__version__)"`.
* The script pins BLAS to one thread per process before importing NumPy (standard
  for MCMC; also avoids a sampler crash seen on some constrained systems).
* Design decisions encoded in the analysis (see Methods): dose enters as
  log10(ppm + 1), centred; each dose-response is quadratic (admits non-monotonic
  shapes); 400 ppm is excluded from the primary fits; SH-SY5Y is excluded from the
  death analysis (harvest artefact) but retained for metabolism; cross-line
  comparison is on within-line responses relative to the 0 ppm control.

## Data

See [`data/DATA_DICTIONARY.md`](data/DATA_DICTIONARY.md) for a full description of
each sheet and its columns. In brief, `Paper_Data_0_3.xlsx` holds the flow
cytometry, resazurin and BCA-protein readouts used by the analysis;
`Paper_Data_Resazurin.xlsx` holds the raw resazurin fluorescence used only by
`--check-resazurin`.

## Citing

If you use this code or data, please cite the manuscript (above) and this
repository — see [`CITATION.cff`](CITATION.cff). On archiving a release to Zenodo
you will receive a DOI; add it to `CITATION.cff` and to the manuscript's data-
availability statement.

## Licence

This repository is dual-licensed:

- **Code** (`co_analysis.py`) — MIT Licence (see [`LICENSE`](LICENSE)).
- **Data** (`data/`) — Creative Commons Attribution 4.0 International, CC BY 4.0
  (see [`data/LICENSE`](data/LICENSE)).

Copyright is held by the University of Sheffield and Sheffield Hallam University
(the authors' institutions). 

## Contact

Alexander G. Fletcher — a.g.fletcher@sheffield.ac.uk
