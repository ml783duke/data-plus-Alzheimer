# ADNI Data Merge & Clean - Workflow Log

## Project Overview
- **Goal**: Merge and clean ADNI multi-source data, anchored on the SOMAscan proteomics dataset
- **Primary Dataset**: `CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx` (737 rows × 7015 columns, ~7000 protein measurements)
- **Working Directory**: `/Users/limai/Desktop/DATA+/ADNI-data/`
- **R Code Directory**: `/Users/limai/Desktop/DATA+/ADNI-data/R-code/`
- **Created**: 2026-06-18

---

## Task Log

| # | Date | Task | Status | Script | Notes |
|---|------|------|--------|--------|-------|
| 0 | 2026-06-18 | Initialize project structure | ✅ Done | `00_reference.R`, `00_reference.md` | Created R-code folder, shared ref files |
| 1 | 2026-06-18 | Survey all 12 data files | ✅ Done | - | Recorded column counts, ID column names, special cases |
| 2 | 2026-06-18 | Examine SOMAscan VISCODE2 distribution | ✅ Done | - | bl=708, NA=7, m24-m72=22. Mostly cross-sectional but not pure. 2 RIDs have duplicate rows |
| 3 | 2026-06-18 | Confirm ID strategy, join type, packages | ✅ Done | - | 4-digit padding, full join, tidyverse |

---

## Key Decisions Made Today

1. **RID standardization**: Pad to 4 digits (e.g., `"69"` → `"0069"`)
2. **Join type**: Full join (keep all rows from both SOMAscan and merge tables)
3. **SOMAscan as anchor**: Merge direction confirmed
4. **R packages**: `readxl`, `dplyr`, `tidyr`, `stringr`

---

## Pending from User
- Which specific variables to merge from each file
- How to handle SOMAscan's 2 duplicate RIDs (keep both rows or filter)
- Whether to keep all timepoints or filter to baseline

---

## Modification History

| Date | What changed | Why |
|------|-------------|-----|
| 2026-06-18 | Initial setup | Project started |
