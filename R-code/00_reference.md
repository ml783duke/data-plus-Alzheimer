# ADNI Data Reference - Key Information

> **Last Updated**: 2026-06-18 (Session 1)  
> **Purpose**: Central record of data rules, conventions, and decisions for this project.  
> Read this FIRST when resuming work after closing the computer.

---

## 1. ID System

### Two Types of IDs

| ID Type | Format | Example | Source |
|---------|--------|---------|--------|
| **PTID** | `XXX_S_XXXX` (site_blinded_subject) | `027_S_4919` | Most demographic/clinical files |
| **RID** | Numeric (4 digits, may lose leading zeros) | `4919` or `69` | Derived from PTID's last 4 digits |

### PTID → RID Conversion Rule
```
RID = last 4 characters of PTID (after the final underscore)
Example: "027_S_4919" → "4919"
         "027_S_0069" → "0069"
```

### ⚠️ Critical Matching Rules
1. **All merges MUST use RID** as the key — never use PTID directly
2. If a file only has PTID, convert to RID first using the rule above, then merge
3. **RID standardization**: ✅ **DECIDED — Pad to 4 digits** (e.g., `"69"` → `"0069"`, `"4919"` → `"4919"`)
   - Always convert RIDs to **character** type (not numeric)
   - Use `stringr::str_pad(rid, width = 4, side = "left", pad = "0")`
4. **⚠️ Merge key = RID + Timepoint** (双重匹配！)
   - 不能只用 RID 匹配——一个病人在纵向数据里可能有多行（不同时间点），只用 RID 会导致错误的全排列
   - 合并时必须同时匹配 **RID** 和 **时间变量**（VISCODE2 / Timepoint）
   - 不同文件的时间列名不同，需要先统一：

### Time Column Mapping (All Files)

| File | Time Column | SOMAscan Equivalent |
|------|------------|-------------------|
| `CruchagaLab_CSF_SOMAscan7k...` (主表) | `VISCODE2` | ← anchor |
| `f10_ptdemog.xlsx` | `VISCODE2` | ✅ 同名，直接匹配 |
| `f06_medhist.xlsx` | `VISCODE2` | ✅ 同名 |
| `f04_npi.xlsx` | `VISCODE2` | ✅ 同名 |
| `f05_npiq.xlsx` | `VISCODE2` | ✅ 同名 |
| `f01_duke_uplc.xlsx` | `VISCODE2` | ✅ 同名 |
| `f03_duke_metabolon.xlsx` | `VISCODE2` | ✅ 同名 |
| `f02a_leiden_hph.xlsx` | `VISCODE2` | ✅ 同名 |
| `f02b_leiden_lph.xlsx` | `VISCODE2` | ✅ 同名 |
| `f07_bhr_baseline.xlsx` | `Timepoint` | ⚠️ 需重命名为 VISCODE2 |
| `f08_bhr_longitudinal.xlsx` | `Timepoint` | ⚠️ 需重命名为 VISCODE2 |
| `f09_bhr_caregiver.xlsx` | `Timepoint` | ⚠️ 需重命名为 VISCODE2 |

> **注意**: BHR 系列文件使用 `Timepoint` 列名而非 `VISCODE2`，合并前需先 `rename(VISCODE2 = Timepoint)` 并检查其值是否与 SOMAscan 的 VISCODE2 值（bl, m24, m36, m48, m60, m72）一致。也可能同时存在日期列（如 `EXAMDATE` / `CollectedDate_DRVD`）可作为辅助校验。

### ⚠️ 防御性编程规则
5. **如果文件没有任何时间列**（VISCODE2 / VISCODE / Timepoint）：**立即报错停下**，打印文件中已有的列名让用户手动检查决定。**绝不自动执行合并**。

### Files with Non-Standard ID Column Names
| File | ID Column Name | Maps To |
|------|---------------|---------|
| `f09_bhr_caregiver.xlsx` | `AboutRID` | RID |
| `f09_bhr_caregiver.xlsx` | `AboutPTID` | PTID |
| `f09_bhr_caregiver.xlsx` | `StudyPartnerID` | Caregiver's ID (separate from subject) |

---

## 2. File Inventory

### Primary Dataset (Anchor)
| New Name | Original File | Rows | Cols | ID Col |
|----------|--------------|------|------|--------|
| *(keep original name)* | `CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx` | 737 | 7015 | `RID` |

> SOMAscan: ~7000 proteins (cols X10000.28 to X9999.1), plus `EXAMDATE`, `VISCODE2`, plate info.
> **VISCODE2**: bl=708, NA=7, m24=5, m36=5, m60=5, m48=4, m72=3. Only 2 RIDs duplicate (1016: m24+m72, 388: m24+m48).

### All Other Data Files

| File | Rows | Cols | ID Cols | Content / Lifestyle Category |
|------|------|------|---------|------|
| `f01_duke_uplc.xlsx` | ? | 122 | `PTID`, `RID` | Duke UPLC metabolites → **Nutrition** |
| `f02a_leiden_hph.xlsx` | ? | 63 | `RID` | Leiden oxylipins HPH → **Nutrition** |
| `f02b_leiden_lph.xlsx` | ? | 67 | `RID` | Leiden oxylipins LPH → **Nutrition** |
| `f03_duke_metabolon.xlsx` | ? | 1359 | `PTID`, `RID` | Duke Metabolon HD4 → **Nutrition**, **Smoking** |
| `f04_npi.xlsx` | ? | 168 | `PTID`, `RID` | Neuropsychiatric Inventory → **Sleep** |
| `f05_npiq.xlsx` | ? | 41 | `PTID`, `RID` | NPI Questionnaire → **Sleep** |
| `f06_medhist.xlsx` | ? | 39 | `PTID`, `RID` | Medical history → **Alcohol**, **Smoking** |
| `f07_bhr_baseline.xlsx` | ? | 72 | `RID`, `PTID` | Baseline health questionnaire → **Social** |
| `f08_bhr_longitudinal.xlsx` | ? | 74 | `RID`, `PTID` | Longitudinal health questionnaire → **Social** |
| `f09_bhr_caregiver.xlsx` | ? | 133 | `AboutRID`, `AboutPTID`, `StudyPartnerID` | Caregiver burden → **Social** |
| `f10_ptdemog.xlsx` | ? | 84 | `PTID`, `RID` | Patient demographics → **Physical activity** |

### Reference / Dictionary
| File | Rows | Cols | Content |
|------|------|------|---------|
| `lifestyle_dictionary.xlsx` | 15 | 6 | **Lifestyle变量字典** — 6大类14个变量：Nutrition(4), Sleep(2), Alcohol(2), Smoking(3), Social(2), Physical activity(1)。列：Category, Variable/Construct, Original file name, New file name, Variable(s) in table, DATADIC meaning |

### Temporary/Lock Files (Ignore)
- `.~CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx`
- `.~adni_paper_aligned_lifestyle_variables.xlsx`

### Reference Dictionaries (`dic_reference/`)
| File | Content |
|------|---------|
| `DATADIC_16Jun2026.xlsx` | ADNI Data Dictionary — all variables across all tables (40,402 rows) |
| `adni_paper_aligned_lifestyle_variables.xlsx` | Original (unfiltered) lifestyle variable mapping from paper |

---

## 3. Key Decisions Log

| # | Decision | Options Considered | Chosen | Date |
|---|----------|-------------------|--------|------|
| 1 | ID standardization strategy | (a) Strip zeros; (b) Pad to 4 digits | **(b) Pad to 4 digits** | 2026-06-18 |
| 2 | Merge join type | (a) Left join (keep only SOMAscan RIDs); (b) Full join (keep all RIDs from both tables) | **(b) Full join** — no data loss from either side | 2026-06-18 |
| 3 | SOMAscan anchor | SOMAscan as the primary anchor table | ✅ Confirmed | 2026-06-18 |
| 4 | R packages | readxl, dplyr, tidyr, stringr | ✅ Confirmed — use standard tidyverse | 2026-06-18 |
| 5 | Timepoints | (a) Baseline only; (b) All timepoints | **PENDING** — user aware SOMAscan has 22 non-bl rows, deciding later | 2026-06-18 |
| 7 | Missing time column | If a file has no VISCODE2/VISCODE/Timepoint, STOP with clear warning and ask user | ✅ Decided — defensive check, not automatic | 2026-06-18 |

---

## 4. R Coding Conventions

- **Language**: R (all scripts)
- **Script naming**: `01_xxx.R`, `02_xxx.R`, etc. (execution order by number)
- **File encoding**: UTF-8
- **Relative paths**: Set working directory to `ADNI-data/` before running
- **Dependencies**: List all `library()` calls at top of each script

---

## 5. Pending Questions for User

1. ~~Which strategy for handling leading zeros in RID?~~ → ✅ Pad to 4 digits
2. Which specific variables from each file should be merged? (User will provide later)
3. ~~Which R packages preferred?~~ → ✅ readxl + dplyr + tidyr + stringr
4. Keep all timepoints or filter to baseline only? (SOMAscan: 708 bl vs 22 non-bl + 7 NA — user reviewing)
5. What to do with SOMAscan duplicate RIDs (RID=1016 m24+m72, RID=388 m24+m48)?
