"""
add_ptau.py — Add pTau-217 (plasma + CSF) and birth year to lifestyle_baseline.csv
====================================================================================

Sources:
  Plasma pTau-217 — Fujirebio/Quanterix (pT217_F) primary; C2N Precivity AD2 (pT217_C2N) fill.
                    Units: pg/mL
  CSF pTau-217    — BSHRI NULISA CNS panel (NPQ, already log-scale). QC-passed only.
  Birth year      — ADNI PTDEMOG (PTDOBYY, year extracted)

New columns:
  ptau217_plasma        — plasma pTau-217, pg/mL
  ptau217_plasma_source — "Fujirebio" or "C2N"
  ptau217_plasma_date   — blood draw date
  age_at_plasma_ptau    — approximate age at plasma draw (collection year − birth year)
  ptau217_csf           — CSF pTau-217, NPQ log-scale (do NOT log-transform again)
  ptau217_csf_date      — CSF collection date
  age_at_csf_ptau       — approximate age at CSF draw (collection year − birth year)
  age_at_ptau           — combined age: plasma date preferred, falls back to CSF date
  birth_year            — year of birth (integer)

Outputs (both in lifestyle-merge-RK/output/):
  lifestyle_baseline_ptau.csv  — full baseline with all new columns (N=4888)
  lifestyle_7var_ptau.csv      — 7 key lifestyle vars + pTau + birth year only

Coverage at baseline (N=4888):
  Fujirebio plasma pTau-217:  1134 (23.2%)  +  C2N fill: ~60 more
  NULISA CSF pTau-217:        1498 (30.6%)
  Birth year:                 4856 (99.3%)
"""

import os
import pandas as pd

LIFESTYLE_OUT_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    "lifestyle-merge-RK", "output"
)
PTAU_DATA_DIR = "/Users/reinakobayashi/Downloads/Tables_08Jul2026"
NULISA_FILE   = "/Users/reinakobayashi/Data+/data/BSHRI_PLA_CSF_NULISA_CNS_16Jun2026.csv"

BASELINE_IN   = os.path.join(LIFESTYLE_OUT_DIR, "lifestyle_baseline.csv")
BASELINE_OUT  = os.path.join(LIFESTYLE_OUT_DIR, "lifestyle_baseline_ptau.csv")
VAR7_OUT      = os.path.join(LIFESTYLE_OUT_DIR, "lifestyle_7var_ptau.csv")

KEY_VARS = ["DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK"]
PTAU_COLS = [
    "ptau217_plasma", "ptau217_plasma_source", "ptau217_plasma_date", "age_at_plasma_ptau",
    "ptau217_csf",    "ptau217_csf_date",      "age_at_csf_ptau",
    "age_at_ptau",
]


def section(title, note=None):
    print("\n" + "─" * 72)
    print(title)
    if note:
        print(note)
    print("─" * 72)


def load_baseline_one_per_rid(fpath, value_col, viscode_val="bl"):
    """
    Load a pTau source file, filter to baseline visit, return one row per RID.
    Keeps the row with the highest non-null value when duplicates exist.
    """
    df = pd.read_csv(fpath)
    df["RID"]     = pd.to_numeric(df["RID"],     errors="coerce")
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    if viscode_val and "VISCODE2" in df.columns:
        df = df[df["VISCODE2"] == viscode_val]

    df = df.sort_values(value_col, ascending=False, na_position="last")
    df = df.drop_duplicates(subset="RID", keep="first")

    keep = ["RID", value_col]
    if "EXAMDATE" in df.columns:
        keep.append("EXAMDATE")
    return df[keep].dropna(subset=["RID"])


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1  Load lifestyle_baseline
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 1: Load lifestyle_baseline.csv")

baseline = pd.read_csv(BASELINE_IN)
baseline["RID"] = baseline["RID"].astype(int)
n_patients = len(baseline)
bl_rids = set(baseline["RID"])
print(f"  {n_patients} patients, {baseline.shape[1]} cols")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2  Plasma pTau-217: Fujirebio (primary) + C2N (fill)
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 2: Plasma pTau-217",
    "Fujirebio pT217_F primary; C2N pT217_C2N fills patients missing Fujirebio."
)

fuji = load_baseline_one_per_rid(
    os.path.join(PTAU_DATA_DIR, "All_Subjects_UPENN_PLASMA_FUJIREBIO_QUANTERIX_08Jul2026.csv"),
    "pT217_F",
).rename(columns={"pT217_F": "ptau217_plasma", "EXAMDATE": "ptau217_plasma_date"})
fuji["ptau217_plasma_source"] = "Fujirebio"

c2n = load_baseline_one_per_rid(
    os.path.join(PTAU_DATA_DIR, "All_Subjects_C2N_PRECIVITYAD2_PLASMA_08Jul2026.csv"),
    "pT217_C2N",
).rename(columns={"pT217_C2N": "ptau217_plasma", "EXAMDATE": "ptau217_plasma_date"})
c2n["ptau217_plasma_source"] = "C2N"

c2n_fill = c2n[~c2n["RID"].isin(set(fuji["RID"])) & c2n["ptau217_plasma"].notna()]
plasma = pd.concat([fuji, c2n_fill], ignore_index=True)

n_fuji  = int((plasma["ptau217_plasma_source"] == "Fujirebio").sum())
n_c2n   = int((plasma["ptau217_plasma_source"] == "C2N").sum())
n_match = plasma[plasma["RID"].isin(bl_rids)]["RID"].nunique()
print(f"  Fujirebio: {n_fuji} | C2N fill: {n_c2n}")
print(f"  Matched to lifestyle_baseline: {n_match}/{n_patients} ({n_match/n_patients*100:.1f}%)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3  CSF pTau-217: NULISA
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 3: CSF pTau-217 — BSHRI NULISA",
    "NPQ units (log-scale, do NOT log-transform again). QC-passed baseline rows only.",
)

nulisa = pd.read_csv(NULISA_FILE, low_memory=False)
nulisa["RID"] = pd.to_numeric(nulisa["RID"], errors="coerce")
nulisa["NPQ"] = pd.to_numeric(nulisa["NPQ"], errors="coerce")

csf_ptau = (
    nulisa[
        (nulisa["Target"] == "pTau-217") &
        (nulisa["SampleMatrixType"] == "CSF") &
        (nulisa["SampleQC"] == "passed") &
        (nulisa["VISCODE2"] == "bl")
    ]
    .dropna(subset=["RID", "NPQ"])
    .sort_values("NPQ", ascending=False)
    .drop_duplicates("RID", keep="first")
    [["RID", "NPQ", "EXAMDATE"]]
    .rename(columns={"NPQ": "ptau217_csf", "EXAMDATE": "ptau217_csf_date"})
)

n_csf_match = csf_ptau[csf_ptau["RID"].isin(bl_rids)]["RID"].nunique()
print(f"  NULISA baseline (QC passed): {len(csf_ptau)} unique patients")
print(f"  Matched to lifestyle_baseline: {n_csf_match}/{n_patients} ({n_csf_match/n_patients*100:.1f}%)")
print(f"  NPQ range: {csf_ptau['ptau217_csf'].min():.2f} – {csf_ptau['ptau217_csf'].max():.2f}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4  Birth year: PTDEMOG
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 4: Birth year — PTDEMOG PTDOBYY")

demog = pd.read_csv(os.path.join(PTAU_DATA_DIR, "All_Subjects_PTDEMOG_08Jul2026.csv"))
demog["RID"] = pd.to_numeric(demog["RID"], errors="coerce")

# PTDOBYY is stored as a date string (e.g. "1944-01-01") — extract year
demog["birth_year"] = pd.to_datetime(demog["PTDOBYY"], errors="coerce").dt.year.astype("Int64")

# One row per patient: PTDOBYY is static so any non-null row is fine
birth = (
    demog.dropna(subset=["RID", "birth_year"])
    .drop_duplicates("RID", keep="first")
    [["RID", "birth_year"]]
)

n_by_match = birth[birth["RID"].isin(bl_rids)]["RID"].nunique()
print(f"  Unique RIDs with birth year: {len(birth)}")
print(f"  Matched to lifestyle_baseline: {n_by_match}/{n_patients} ({n_by_match/n_patients*100:.1f}%)")
print(f"  Year range: {int(birth['birth_year'].min())} – {int(birth['birth_year'].max())}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5  Merge all → lifestyle_baseline_ptau.csv
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 5: Merge all → lifestyle_baseline_ptau.csv")

result = (
    baseline
    .merge(plasma[["RID", "ptau217_plasma", "ptau217_plasma_source", "ptau217_plasma_date"]], on="RID", how="left")
    .merge(csf_ptau[["RID", "ptau217_csf", "ptau217_csf_date"]], on="RID", how="left")
    .merge(birth, on="RID", how="left")
)

# Age at collection = collection year − birth year (approximate; only birth year available)
_plasma_year = pd.to_datetime(result["ptau217_plasma_date"], errors="coerce").dt.year
_csf_year    = pd.to_datetime(result["ptau217_csf_date"],    errors="coerce").dt.year

result["age_at_plasma_ptau"] = _plasma_year - result["birth_year"]
result["age_at_csf_ptau"]    = _csf_year    - result["birth_year"]
# Combined: plasma preferred; falls back to CSF when plasma date is missing
result["age_at_ptau"]        = _plasma_year.fillna(_csf_year) - result["birth_year"]

n_plasma = result["ptau217_plasma"].notna().sum()
n_csf    = result["ptau217_csf"].notna().sum()
n_both   = (result["ptau217_plasma"].notna() & result["ptau217_csf"].notna()).sum()
n_by     = result["birth_year"].notna().sum()

print(f"  Output: {len(result)} rows × {result.shape[1]} cols")
print(f"  ptau217_plasma:  {n_plasma}/{n_patients} ({n_plasma/n_patients*100:.1f}%)")
print(f"  ptau217_csf:     {n_csf}/{n_patients} ({n_csf/n_patients*100:.1f}%)")
print(f"  Both pTau:       {n_both}/{n_patients} ({n_both/n_patients*100:.1f}%)")
print(f"  birth_year:      {n_by}/{n_patients} ({n_by/n_patients*100:.1f}%)")

result.to_csv(BASELINE_OUT, index=False)
print(f"  Saved: {BASELINE_OUT}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6  7-variable filtered version → lifestyle_7var_ptau.csv
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 6: 7-variable subset → lifestyle_7var_ptau.csv",
    f"Key vars: {KEY_VARS}"
)

present = [v for v in KEY_VARS if v in result.columns]
missing = [v for v in KEY_VARS if v not in result.columns]
if missing:
    print(f"  WARNING: not found in result: {missing}")

keep_cols = ["RID"] + present + PTAU_COLS + ["birth_year"]
var7 = result[keep_cols].copy()

# Drop rows where ALL lifestyle vars AND all pTau values are missing
any_data_cols = present + ["ptau217_plasma", "ptau217_csf"]
var7 = var7.dropna(subset=any_data_cols, how="all")

print(f"  Rows with at least one value: {len(var7)}")
print(f"  Coverage per column:")
for col in keep_cols[1:]:
    n = var7[col].notna().sum()
    print(f"    {col}: {n}/{len(var7)} ({n/len(var7)*100:.1f}%)")

var7.to_csv(VAR7_OUT, index=False)
print(f"  Saved: {VAR7_OUT}")

print("\n" + "═" * 72)
print("DONE — add_ptau.py")
print(f"  {BASELINE_OUT}")
print(f"    → {result.shape[1]} cols ({len(result)} patients)")
print(f"  {VAR7_OUT}")
print(f"    → {len(keep_cols)} cols ({len(var7)} patients with ≥1 value)")
print("═" * 72)
