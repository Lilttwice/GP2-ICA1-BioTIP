# GP2 ICA1 — BioTIP analysis (GSE141259 / GSE87038)

GP2 Data Analysis ICA: BioTIP case study and lung injury scRNA-seq analyses for ZJE Y2 BMI sutdents.

## Setup

1. Open `Your_ICA.Rproj` in RStudio.
3. Install packages: BioTIP, Seurat, GEOquery, MouseGastrulationData (for case study).
4. Run `source("STEP1_check_environment.R")` to verify the environment.

## Data

Raw GEO files are **not** in this repo. Download locally:

```r
source("CODE/M030_gse141259_preprocessing/03_code/download_gse141259_data.R")
source("CODE/M030_gse141259_preprocessing/03_code/preprocess_gse141259.R")
```


## Main scripts (project root)

| Script | Purpose |
|--------|---------|
| `STEP1_check_environment.R` | Environment check |
| `STEP4_install_and_run_case_study.R` | GSE87038 BioTIP case study |
| `STEP6_run_lung_AB_figures.R` | Lung Analyses A & B + figures |
| `STEP6B_fix_analysisB_figures.R` | Fix Analysis B outputs |

Figures for the report: `figures_for_word/`

## Module layout (`CODE/`)

| Folder | Purpose |
|--------|---------|
| `M010_biotip_case_study` | GSE87038 validation |
| `M030_gse141259_preprocessing` | Download + preprocess lung data |
| `M040_analysisA_timepoint_biotip` | Analysis A |
| `M060_analysisB_cellstate_biotip` | Analysis B |
| `M070_compare_and_comment` | Compare A/B |
| `M080_report_assembly` | Report figures/tables |

## Report

- Draft structure: `Your_report_skeleton.md`
- Manuscript text: `GP2_ICA1_终稿正文.md`
