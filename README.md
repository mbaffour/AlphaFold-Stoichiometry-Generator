# AlphaFold Stoichiometry Generator

A Shiny app for generating AlphaFold multimer input files with systematic stoichiometry exploration.

## What It Does

- **Input** protein sequences via paste or FASTA upload (any number of proteins)
- **Set stoichiometry ranges** per protein (min/max copy numbers)
- **Generate combinations** — all combos, random subsets, or homomer series
- **Preview** AlphaFold3 JSON, multimer FASTA, and text summaries
- **Download** individual models or bulk ZIP with manifest

## Install & Run

```r
# Install dependencies (one time)
install.packages(c("shiny", "DT", "jsonlite"))

# Run the app
shiny::runApp("app.R")
```

No Bioconductor dependencies. No external APIs. Three CRAN packages. Works offline.

## Usage

1. **Paste or upload** your protein sequences in FASTA format
2. **Click Parse** to load them
3. **Set min/max copy numbers** for each protein in the sidebar
4. **Choose generation mode:**
   - *All combinations* — full cartesian product of copy number ranges
   - *Random subset* — sample N random stoichiometries
   - *Quick homomer series* — each protein independently (monomer → multimer)
5. **Click Generate Models**
6. **Preview** any model's JSON/FASTA in the Preview tab
7. **Download** from the Download tab

## Output Formats

### AlphaFold3 JSON
Ready to paste into the AlphaFold3 server. Uses the `proteinChain` format with `count` field.

### Multimer FASTA
Each chain as a separate entry (e.g., `>ProteinA_copy1`, `>ProteinA_copy2`).

### Manifest CSV
Tracks all generated models with chain counts, residue totals, and warnings.

## Residue Limit Warning

Models exceeding 3,000 total residues are flagged — AlphaFold performance degrades significantly above this threshold.

## Deploy

For lab-wide access, deploy to:
- **shinyapps.io** (free tier works for small labs)
- **Texas A&M Shiny Server** (if available)
- **Any R server** with `shiny::runApp()`
