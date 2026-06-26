# Lifestyle & Protein Associations with PTAU in ADNI

## Research Progress Report — 2026-06-25

---

## 1. Research Question & Framework

### Core Question

> How do **lifestyle factors** and **CSF proteins** relate to **PTAU** (a core Alzheimer's pathology marker)?

### Conceptual Model

```
┌──────────────────────────────────────────────────────┐
│                    LIFESTYLE                          │
│  Nutrition · Smoking · Alcohol · Sleep · Social      │
└────────┬──────────────────────────┬──────────────────┘
         │                          │
   Phase 2 (next)             Phase 1 (done)
   Lifestyle → Protein         Lifestyle → PTAU
         │                     Result: NO direct effect
         ↓                          │
┌────────────────┐                   │
│    PROTEIN      │←─────────────────┘
│   (7000+)       │
└───────┬─────────┘
        │
  Phase 0 (done)
  Protein → PTAU
  Result: 485 significant
        │
        ↓
┌───────────┐
│   PTAU    │
└───────────┘
```

### Four-Phase Research Plan

| Phase  | Question                                         | Method                          |   Status   |
| ------ | ------------------------------------------------ | ------------------------------- | :--------: |
| **0a** | Which proteins correlate with PTAU?              | Spearman correlation + FDR      |  ✅ Done   |
| **0b** | Which proteins' age-trajectory tracks PTAU?      | LOESS trend correlation         |  ✅ Done   |
| **0c** | Does the AD spectrum change these relationships? | LOESS trend (AD only)           |  ✅ Done   |
| **1**  | Does lifestyle directly relate to PTAU?          | Linear regression + AGE/SEX     |  ✅ Done   |
| **2**  | Which proteins are lifestyle-sensitive?          | Lifestyle × 7000 protein screen |  🔲 Next   |
| **3**  | Does lifestyle affect PTAU through proteins?     | Mediation analysis              | 🔲 Planned |
| **4**  | Do proteins modify lifestyle-PTAU links?         | Interaction analysis            | 🔲 Planned |

---

## 2. Data Overview

| Item            | Detail                                                                                                                                              |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Cohort**      | ADNI baseline (VISCODE2 = "bl")                                                                                                                     |
| **Sample size** | 495 subjects (baseline + PTAU + AGE); AD spectrum subset: 368                                                                                        |
| **Proteomics**  | SomaScan 7k: 7,005 plasma proteins (RFU, log2-transformed)                                                                                          |
| **Lifestyle**   | 43 variables: nutrition (DHA/EPA/fatty acids/oxylipins), smoking (status/pack-years/cotinine), alcohol, sleep (NPI/NPIQ), social, physical activity |
| **AD Spectrum** | CN=127, EMCI=125, LMCI=146, AD=97                                                                                                                   |
| **Main file**   | `master_data.xlsx`: 737 rows × 7,072 columns                                                                                                        |

---

## 3. Completed Analyses & Key Results

### 3.1 Protein ↔ PTAU: Cross-Sectional Correlation

**Method**: Spearman rank correlation, log2(protein) vs log2(PTAU), Benjamini-Hochberg FDR correction

**Results**:

- **7,005 proteins tested** (3 removed for >20% missing)
- **485 proteins significant at FDR < 0.05** (6.9%)
- **rho range**: -0.76 to +0.86

**Top positive correlates (highest rho)**:

| Gene   | Protein                      | rho   | FDR      |
| ------ | ---------------------------- | ----- | -------- |
| YWHAG  | 14-3-3 protein gamma         | 0.861 | 2.6e-143 |
| YWHAE  | 14-3-3 protein epsilon       | 0.839 | 1.5e-128 |
| PPP3R1 | Calcineurin subunit B type 1 | 0.839 | 4.6e-128 |
| GAP43  | Neuromodulin                 | 0.821 | 2.8e-117 |
| YWHAB  | 14-3-3 protein beta/alpha    | 0.818 | 4.3e-115 |

> **Key finding**: The **14-3-3 protein family** (YWHAG, YWHAE, YWHAB, YWHAH) dominates the top PTAU correlates. 14-3-3 proteins are established CSF biomarkers for neurodegeneration.

---

### 3.2 Protein ↔ PTAU: Age-Trend Tracking (Full Population)

**Method**: LOESS smoothing of log2(PTAU) ~ AGE and log2(protein) ~ AGE; Spearman correlation between the two fitted curves = "trend similarity score." 495 baseline subjects.

**Key context**: PTAU itself **strongly increases with age** (Spearman rho = 0.797, p = 7.3e-110).

**Top 10 proteins whose age-trajectory best matches PTAU**:

| Rank | Gene              | Protein                                        | Trend ρ | MAE_z |
| :--: | ----------------- | ---------------------------------------------- | ------- | ----- |
|  1   | **OMP**           | Olfactory marker protein                       | 0.990   | 0.074 |
|  2   | **NCKIPSD**       | NCK-interacting protein with SH3 domain        | 0.989   | 0.184 |
|  3   | **LIN7B**         | Protein lin-7 homolog B                        | 0.988   | 0.097 |
|  4   | **CNRIP1**        | CB1 cannabinoid receptor-interacting protein 1 | 0.983   | 0.092 |
|  5   | **CRISP2**        | Cysteine-rich secretory protein 2              | 0.978   | 0.160 |
|  6   | **UBE2Q2**        | Ubiquitin-conjugating enzyme E2 Q2             | 0.975   | 0.194 |
|  7   | **PPP3CA/PPP3R1** | Calcineurin                                    | 0.975   | 0.313 |
|  8   | **GLIPR2**        | GAPR-1                                         | 0.974   | 0.258 |
|  9   | **STOML2**        | Stomatin-like protein 2                        | 0.972   | 0.152 |
|  10  | **IL18RAP**       | IL-18 Rb                                       | 0.972   | 0.113 |

> **Key finding**: These are NOT the same proteins found by cross-sectional correlation. Only **2 out of 60 overlap** between the two approaches — demonstrating that "age-trajectory tracking" and "cross-sectional correlation" capture fundamentally different biology.

---

### 3.3 Protein ↔ PTAU: AD Spectrum Only (EMCI+LMCI+AD)

**Method**: Same LOESS trend analysis, restricted to 368 cognitively impaired subjects, excluding 127 CN.

**Top 10 proteins in the AD spectrum**:

| Rank | Gene       | Protein | Trend ρ | Rank in Full Pop |
| :--: | ---------- | ------- | ------- | :--------------: |
|  1   | **CNOT9**  | RCD1    | 0.989   |      #256 ⬆      |
|  2   | **LIN7B**  | LIN7B   | 0.987   |        #3        |
|  3   | **SH3GL2** | SH3G2   | 0.987   |      #149 ⬆      |
|  4   | **YWHAE**  | 14-3-3E | 0.986   |      #13 ⬆       |
|  5   | **FAM20A** | FA20A   | 0.983   |      #309 ⬆      |
|  6   | **PTPN11** | SHP-2   | 0.982   |      #141 ⬆      |
|  7   | **TMX1**   | TXND1   | 0.980   |      #251 ⬆      |
|  8   | **AARS1**  | SYAC    | 0.980   |      #125 ⬆      |
|  9   | **CNRIP1** | CB032   | 0.979   |        #4        |
|  10  | **PSMG3**  | PSMG3   | 0.978   |      #276 ⬆      |

> **Key findings**:
>
> - Only **3 / 30 proteins overlap** between full-population and AD-spectrum top lists
> - **OMP** (Full #1) drops to **#475** in AD → its PTAU-tracking relationship may be disrupted by AD pathology
> - **Calcineurin** (PPP3CA/PPP3R1, Full #7) plummets to **#1,958** in AD → phosphatase dysfunction in AD
> - **CNOT9, SH3GL2** emerge as new top trackers in the AD continuum (ranked #256 and #149 in full pop)
> - **YWHAE (14-3-3 epsilon)** is the only protein consistently top-ranked across ALL three analyses (correlation #2, full trend #13, AD trend #4)

---

### 3.4 Lifestyle → PTAU: Direct Associations

**Method**: `log2(PTAU) ~ lifestyle_var + AGE + SEX`, 25 lifestyle variables tested, FDR corrected.

**Results**: **No lifestyle variable reached FDR significance** (all FDR > 0.19).

Top raw associations (none survive multiple-testing correction):

| Variable           | β (per 1-SD) | Raw P | FDR  |
| ------------------ | :----------: | :---: | :--: |
| DHA oxylipin (LPH) |    -0.085    | 0.020 | 0.20 |
| EPA oxylipin (HPH) |    -0.083    | 0.023 | 0.20 |
| Cotinine           |    -0.071    | 0.024 | 0.20 |
| Hydroxycotinine 2  |    -0.067    | 0.032 | 0.20 |

> **Key finding**: Lifestyle factors do **not** show direct, independent associations with PTAU after adjusting for age and sex. This suggests that any lifestyle effect on PTAU is likely **indirect** — mediated through changes in protein expression or metabolism. This motivates Phase 2–3.

---

## 4. Cross-Analysis Consistency: The 14-3-3 Story

**YWHAE (14-3-3 epsilon)** is the standout protein across all completed analyses:

| Analysis                    | YWHAE Performance                         |
| --------------------------- | ----------------------------------------- |
| Cross-sectional correlation | **#2** (rho = 0.839, FDR = 1.5e-128)      |
| Full-population age trend   | **#13** (trend ρ = 0.969)                 |
| AD-spectrum age trend       | **#4** (trend ρ = 0.986) ⬆ improves in AD |
| PTAU-associated Excel list  | Listed as **Positive**                    |

Other 14-3-3 family members (YWHAG, YWHAB, YWHAH) also rank highly across analyses, confirming the 14-3-3 family as a **core PTAU-associated protein module**.

---

## 5. Summary of Completed vs. Planned Work

```
✅ DONE:
   ├── Phase 0a: Protein-PTAU cross-sectional correlation (485 FDR<0.05)
   ├── Phase 0b: Protein-PTAU age-trend tracking, full pop (7005 ranked)
   ├── Phase 0c: Protein-PTAU age-trend tracking, AD spectrum (7003 ranked)
   └── Phase 1:  Lifestyle → PTAU direct (no significant findings)

🔲 NEXT:
   ├── Phase 2:  Lifestyle → Protein screen
   │             For each significant lifestyle var × 7000 proteins:
   │             log2(protein) ~ lifestyle + AGE + SEX, FDR corrected
   │             → Identify "lifestyle-sensitive proteins"
   │
   ├── Phase 3:  Mediation: Lifestyle → Protein → PTAU
   │             (on Phase 1+2 significant combos only)
   │
   └── Phase 4:  Interaction: Lifestyle × Protein → PTAU
```

---

## 6. Output Files

All results are in `output/`:

| Directory                    | Content                                          |
| ---------------------------- | ------------------------------------------------ |
| `ptau_protein_correlation/`  | Volcano plot, top30 bar plots, CSVs              |
| `ptau_protein_age_trend/`    | Combined trend plot (top10 + PTAU), ranking CSVs |
| `ptau_protein_age_trend_AD/` | Same for AD spectrum only                        |
| `lifestyle_ptau/`            | Forest plot, association CSVs                    |

Scripts: `R_code/01–04_*.R`
