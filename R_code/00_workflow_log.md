# ADNI Multi-Omics + Lifestyle Research — Workflow Log

## Project Overview

- **Goal**: Understand how lifestyle factors and plasma proteins relate to PTAU (Alzheimer's pathology marker) in aging and AD
- **Data**: ADNI baseline samples, SomaScan 7k proteomics + 14 lifestyle variables + clinical/demographic
- **Main data file**: `master_data.xlsx` (737 rows × 7029 cols)
- **Working directory**: `/Users/limai/Desktop/DATA+/ADNI-data/`
- **Created**: 2026-06-18 | **Last updated**: 2026-06-25

---

## Project Structure

```
DATA+/ADNI-data/
├── master_data.xlsx              # Main merged dataset
├── R_code/                       # All R scripts
│   ├── 00_reference.R            # Shared utilities, file paths, ID functions
│   ├── 00_reference.md           # Reference documentation
│   ├── 00_workflow_log.md        # This file
│   ├── 01_ptau_protein_correlation.R   # Protein-PTAU Spearman (cross-sectional)
│   ├── 02_ptau_protein_age_trend.R     # Protein-PTAU age-trend tracking (full pop)
│   ├── 03_ptau_protein_age_trend_AD.R  # Protein-PTAU age-trend tracking (AD spectrum)
│   ├── 04_lifestyle_ptau.R             # [PLANNED] Lifestyle → PTAU direct
│   ├── 05_lifestyle_protein.R          # [PLANNED] Lifestyle → Protein screening
│   ├── 06_mediation.R                  # [PLANNED] Lifestyle → Protein → PTAU mediation
│   ├── 07_interaction.R                # [PLANNED] Lifestyle × Protein → PTAU
│   ├── add_new_ptau_compare.R          # Add NULISA pTau-217 to master
│   ├── merge 2 somascan to master.R    # Merge SomaScan + clinical into master
│   └── check if 2 somascan is the same (it is).R  # Validation script
├── protein_raw_data/             # Proteomics data
│   ├── CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx
│   ├── somascan with gene ID.txt
│   ├── protein dict.csv          # SeqId → gene symbol mapping
│   └── BSHRI_PLA_CSF_NULISA_CNS_22Jun2026.xlsx  # NULISA pTau-217
├── lifestyle_raw_data/           # Lifestyle data (f01-f10)
│   ├── lifestyle_dictionary.xlsx # 14 lifestyle variables catalog
│   ├── f01_duke_uplc.xlsx       # Nutrition: DHA, EPA, homocysteine
│   ├── f02a_leiden_hph.xlsx     # Nutrition: oxylipins (high confidence)
│   ├── f02b_leiden_lph.xlsx     # Nutrition: oxylipins (low confidence)
│   ├── f03_duke_metabolon.xlsx  # Nutrition: fatty acids, cotinine
│   ├── f04_npi.xlsx             # Sleep: NPI sleep items
│   ├── f05_npiq.xlsx            # Sleep: NPIQ sleep items
│   ├── f06_medhist.xlsx         # Smoking + Alcohol history
│   ├── f07_bhr_baseline.xlsx    # Social: reduced activities
│   ├── f08_bhr_longitudinal.xlsx # Social: reduced activities (longitudinal)
│   ├── f09_bhr_caregiver.xlsx   # Social: loneliness
│   └── f10_ptdemog.xlsx         # Physical activity: occupation
├── dic_reference/                # Data dictionaries
└── output/                       # All analysis outputs
    ├── ptau_protein_correlation/      # Output of script 01
    ├── ptau_protein_age_trend/        # Output of script 02
    └── ptau_protein_age_trend_AD/     # Output of script 03
```

---

## Completed Work

### Phase 0: Data Preparation

| # | Date | Task | Script | Key Results |
|---|------|------|--------|-------------|
| 0.1 | 2026-06-18 | Project init, survey all data files | `00_reference.R` | RID strategy (4-digit pad), full join, SomaScan as anchor |
| 0.2 | 2026-06-18 | Verify two SomaScan sources are identical | `check if 2 somascan is the same (it is).R` | Confirmed: long-format txt = wide-format Excel (log2 values) |
| 0.3 | 2026-06-18 | Merge SomaScan + clinical variables → master_data.xlsx | `merge 2 somascan to master.R` | 737 × 7029, all relative paths |
| 0.4 | 2026-06-20 | Restructure project folders | — | Created `protein_raw_data/`, `lifestyle_raw_data/`, moved files |
| 0.5 | 2026-06-20 | Fix all R file paths to relative | — | All 5 scripts use relative paths from project root |
| 0.6 | 2026-06-20 | Restructure output directory | — | `output/` at project root, subfolders per analysis |
| 0.7 | 2026-06-20 | Add NULISA pTau-217 to master | `add_new_ptau_compare.R` | 495 baseline with PTAU, Spearman PTAU vs new_ptau rho |
| 0.8 | 2026-06-25 | Fix gene annotation bug in 01 + 02 | — | Merge key: `Analytes` column (not `SeqId`) — has `X` prefix |

### Phase 1: Protein ↔ PTAU Discovery

| # | Date | Task | Script | Key Results |
|---|------|------|--------|-------------|
| 1.1 | 2026-06-20 | PTAU-protein Spearman correlation | `01_ptau_protein_correlation.R` | 7005 proteins tested, 485 FDR<0.05, volcano + bar plots, top30 pos/neg with gene names |
| 1.2 | 2026-06-25 | PTAU-protein age-trend tracking (full pop) | `02_ptau_protein_age_trend.R` | 495 subjects. LOESS trend correlation. Top: OMP, NCKIPSD, LIN7B, CNRIP1, CRISP2. Only 2/60 overlap with 1.1's FDR-significant list |
| 1.3 | 2026-06-25 | PTAU-protein age-trend tracking (AD spectrum) | `03_ptau_protein_age_trend_AD.R` | 368 subjects (EMCI+LMCI+AD). Top: CNOT9, LIN7B, SH3GL2, YWHAE, FAM20A. Overlap with full pop top30: only 3 proteins. 14-3-3E (YWHAE) consistently top-ranked across all analyses |

### Key Findings Summary

- **PTAU increases with age**: Spearman rho = 0.80 (full pop) / 0.80 (AD spectrum), both p ≈ 0
- **Cross-sectional correlation vs age-trend tracking are fundamentally different questions**: Only 2/60 proteins overlap between the two approaches
- **YWHAE (14-3-3 epsilon)** is the only protein that performs well in ALL three analyses (correlation, full-pop trend, AD trend)
- **AD spectrum reshuffles the protein-PTAU age relationship**: Full pop #1 (OMP) drops to #475 in AD; new leaders emerge (CNOT9 #1, SH3GL2 #3)
- **14-3-3 family** (YWHAE, YWHAG, YWHAB, YWHAH) strongly associated with PTAU across multiple frameworks

---

## Research Framework: Lifestyle + Protein → PTAU

### Conceptual Model

```
                    ┌──────────────┐
                    │   LIFESTYLE   │
                    │  (14 vars)    │
                    │ Nutrition     │
                    │ Sleep         │
                    │ Smoking       │
                    │ Alcohol       │
                    │ Social        │
                    │ Phys Activity │
                    └───┬──────┬────┘
                        │      │
              Phase 2   │      │  Phase 1
              (screen)  │      │  (direct)
                        ↓      ↓
                    ┌──────────────┐
                    │   PROTEIN     │←──────┐
                    │  (7000+)      │       │ Phase 3 (mediation)
                    └───────┬───────┘       │
                            │               │
                 Phase 1    │               │
                 (done)     ↓               │
                        ┌──────────┐        │
                        │   PTAU    │←──────┘
                        └──────────┘
                                 ↑
              Phase 4 (interaction): Lifestyle × Protein → PTAU
```

### Planned Analyses

| Phase | Script | Question | Method | Input | Output |
|-------|--------|----------|--------|-------|--------|
| **1** | `04_lifestyle_ptau.R` | Which lifestyle factors directly relate to PTAU? | `log2(PTAU) ~ lifestyle + AGE + SEX` for each of 14 vars | master + lifestyle | Table of 14 associations |
| **2** | `05_lifestyle_protein.R` | Which proteins are affected by lifestyle? | `log2(prot) ~ lifestyle + AGE + SEX` for significant lifestyle vars × 7000 proteins, FDR corrected | master + lifestyle | Heatmap: lifestyle × protein sensitivity |
| **3** | `06_mediation.R` | Does lifestyle affect PTAU through proteins? | Mediation: Lifestyle → Protein → PTAU (Baron-Kenny, on phase 1+2 significant combos only) | master + lifestyle | Mediation proportions, key pathways |
| **4** | `07_interaction.R` | Do proteins modify lifestyle-PTAU relationship? | `log2(PTAU) ~ lifestyle × log2(prot) + AGE + SEX` | master + lifestyle | Significant interactions (protective / synergistic) |

### Scope Narrowing Strategy

```
Phase 1: 14 lifestyle vars → ~5 significant
Phase 2: 5 lifestyle × 7000 proteins → ~200 lifestyle-sensitive proteins (FDR<0.05)
Phase 3: 5 × 200 → ~10 significant mediation pathways
Phase 4: 5 × 200 → ~5 significant interactions
```

Each phase reduces the search space for the next, preventing multiple-testing overload and keeping results interpretable.

### Lifestyle Variables Catalog

| Category | Variable | Source | Columns |
|----------|----------|--------|---------|
| Nutrition | n-3 PUFA (DHA, EPA) | f01_duke_uplc | DHA, EPA |
| Nutrition | Homocysteine (Hcy) | f01_duke_uplc | HCys |
| Nutrition | Oxylipins (EPA/DHA derived) | f02a/b_leiden | BSH_FA20.5_w3, BSH_FA22.6_w3 |
| Nutrition | Fatty acid biomarkers | f03_duke_metabolon | X424, X439, X10000893… |
| Sleep | NPI sleep domain | f04_npi | NPIK (sleep items) |
| Sleep | NPIQ sleep items | f05_npiq | NPIK, NPIKSEV |
| Smoking | Smoking status | f06_medhist | MH16SMOK |
| Smoking | Pack-years | f06_medhist | MH16ASMOK, MH16BSMOK… |
| Smoking | Cotinine (biomarker) | f03_duke_metabolon | X848, X100002717… |
| Alcohol | Abuse history | f06_medhist | MH14ALCH |
| Alcohol | Drinks/day | f06_medhist | MH14AALCH |
| Social | Reduced social activities | f07_bhr_baseline, f08 | QID46_9, QID50_9 |
| Social | Loneliness | f09_bhr_caregiver | QID44_1_3, QID44_2_3 |
| Phys Activity | Occupation (proxy) | f10_ptdemog | PTWORK, PTWORKHS |

---

## Key Decisions

1. **RID standardization**: Pad to 4 digits (e.g., `"69"` → `"0069"`)
2. **Join type**: Full join (keep all rows from both sides)
3. **SomaScan as anchor**: Merge direction confirmed
4. **Baseline only** for all protein analyses (`VISCODE2 == "bl"`)
5. **Log2 transformation** for PTAU and all protein RFU values
6. **>20% missing filter** for proteins (within analysis set)
7. **All file paths relative** from project root for collaborator compatibility
8. **Gene annotation**: Use `Analytes` column from protein dict (format: `X10000.28`) as merge key
9. **LOESS span**: 0.75 default, auto-widen for small N (`max(0.75, 15/n_valid)`)
10. **AD spectrum definition**: DX in {EMCI, LMCI, AD}, excluding CN

---

## Modification History

| Date | What changed | Why |
|------|-------------|-----|
| 2026-06-18 | Initial setup | Project started |
| 2026-06-20 | Folder restructure, path fixes, new_ptau added | Collaborator compatibility |
| 2026-06-25 | Complete workflow rewrite; added phases 1-3 (done) + phases 4-7 (planned) | Research framework clarified |
