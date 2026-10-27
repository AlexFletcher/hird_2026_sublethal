# Data dictionary

Two workbooks accompany the analysis. `Paper_Data_0_3.xlsx` is the consolidated
dataset read by the analysis; `Paper_Data_Resazurin.xlsx` holds the raw resazurin
fluorescence and is used only by the `--check-resazurin` provenance check.

The sheets were produced and laid out by the experimenters for human reading, so
some are "wide" (cell lines side by side, replicates in blocks). The loaders in
`co_analysis.py` parse these layouts into tidy tables; the descriptions below match
what those loaders expect.

## Cell lines

| Code | Cell line | Origin | Status |
|------|-----------|--------|--------|
| U87 | U87-MG | Glioblastoma (astrocytic) | Transformed |
| SH-SY5Y | SH-SY5Y | Neuroblastoma (neuronal), undifferentiated | Transformed |
| HA | Human astrocytes | Primary astrocyte | Primary |
| HUVEC | Human umbilical vein endothelial cells | Primary endothelial | Primary |

CO dose ladder (nominal gas-phase): 0, 10, 20, 40, 100, 120, 400 ppm; 24 h exposure.
400 ppm is excluded from the primary model fits.

---

## `Paper_Data_0_3.xlsx`

### Sheet `(1-2) Flow Data`  — *used by the analysis*
Annexin V / propidium iodide (PI) flow cytometry. Cell lines are laid out side by
side; each biological replicate is a band labelled `<LINE> One/Two/Three`. Under
each band label, a header row then seven dose rows. Columns per block:

| Column | Meaning |
|--------|---------|
| (dose) | CO concentration (ppm) for the row |
| Cells | total cells in the live/dead gate |
| Live | Annexin V− / PI− (viable; reference category) |
| Dead | dead count from the live/dead gate |
| Annexin V | Annexin V+ / PI− (early apoptotic) |
| Double + | Annexin V+ / PI+ (late apoptotic) |
| PI | Annexin V− / PI+ (necrotic) |

A "Supplementary Info" controls block sits to the right of each band (single-stain
compensation controls: Negative / AV / PI; and staurosporine positive controls).
The analysis maps the four quadrants to live/early/late/necrotic, sums them as the
multinomial total, and computes quadrant fractions.

### Sheet `(3) Resazurin`  — *used by the analysis*
Processed (background-corrected) resazurin fluorescence intensity (FI), the
metabolic-activity readout. Per cell line, seven dose rows; three runs (biological
replicates, different passages) of six wells each. Values are corrected FI
(raw − mean media − mean blank). Non-positive / flagged wells are dropped by the
loader; within each (line, run) the analysis expresses activity relative to the
0 ppm control.

### Sheet `(3) Resazurin Raw`
Raw resazurin FI underlying `(3) Resazurin`, in the consolidated workbook's layout.
(The analysis reads the processed sheet; the standalone raw workbook below is what
the provenance check uses.)

### Sheets `Protein Concentration Cancer` / `Protein Concentration Healthy`  — *used by the analysis*
Back-calculated total protein concentration (µg/mL) by BCA assay — the biomass
normaliser. *Cancer* holds U87 and SH-SY5Y; *Healthy* holds HA and HUVEC. Per line,
seven dose rows; six 1-in-10 readings per row (3 runs × 2 technical reads), followed
by a 1-in-5 dilution block used as a QC cross-check. Absolute levels differ between
lines (biological variability); the analysis works within-line, relative to the
0 ppm control.

### Sheet `(5) Cell Meta-data`
Per cell line: supplier, catalogue number, passage range (source for manuscript
Table 1).

### Sheets `(4) Design` / `(4) Checklist`
Experimental design summary and a data-collection checklist (documentation; not read
by the analysis).

---

## `Paper_Data_Resazurin.xlsx`  — *used only by `--check-resazurin`*

Raw resazurin fluorescence in the original Healthy/Cancer-sheet layout.

- Sheet `Healthy` — HA, HUVEC
- Sheet `Cancer` — U87, SH-SY5Y

Each sheet contains one block per (dose, line). A block header reads
`<dose>ppm - <line>`, followed by rows:

| Row | Meaning |
|-----|---------|
| Res | raw resazurin FI for each well |
| Media | media-only wells (cellular autofluorescence) |
| Blanks | resazurin-in-PBS wells (reagent background) |
| Corrected R | the stored corrected value (= Res − mean Media − mean Blank) |

Three runs of six wells each. `--check-resazurin` recomputes the correction from
Res/Media/Blanks and verifies it against the stored "Corrected R" row. Wells flagged
with a red fill (e.g. a contaminated well) are treated as excluded.
