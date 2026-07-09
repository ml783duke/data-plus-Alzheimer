"""
add_dx_apoe.py — Add cognitive diagnosis and APOE genotype
===========================================================

Reads lifestyle_baseline_ptau.csv and lifestyle_7var_ptau.csv,
appends three new columns, and overwrites the same files.

New columns:
  dx_entry        — diagnosis at study entry (from ADNI Study Entry file)
                    Values: CN, MCI, EMCI, LMCI, SMC, AD, Patient
  dx_bl           — DXSUM diagnosis at baseline (or latest available visit)
                    Values: CU, MCI, Dementia  (null if no DXSUM row found)
  dx_bl_visit     — which DXSUM visit was used ("bl" or actual viscode if fallback)
  apoe_genotype   — APOE genotype from APOERES (e.g. "3/3", "3/4", "4/4")

Sources:
  Study entry:  Study_Entry_09Jul2026.csv
  DXSUM:        DXSUM_08Jul2026.csv
  APOERES:      APOERES_08Jul2026.csv

Run after add_ptau.py:
  python3 add_dx_apoe.py
"""

import os
import pandas as pd

# Specify data path
LIFESTYLE_OUT_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "output")

BASELINE_PTAU = os.path.join(LIFESTYLE_OUT_DIR, "lifestyle_baseline_ptau.csv")
VAR7_PTAU = os.path.join(LIFESTYLE_OUT_DIR, "lifestyle_7var_ptau.csv")

STUDY_ENTRY_FILE = ""
DXSUM_FILE = ""
APOERES_FILE = ""

DX_MAP = {1: "CU", 2: "MCI", 3: "Dementia"}

NEW_COLS = ["dx_entry", "dx_bl", "dx_bl_visit", "apoe_genotype"]

SEP = "─" * 72


def section(title):
    print(f"\n{SEP}\n{title}\n{SEP}")


def pct(n, total):
    return f"{n}/{total} ({n/total*100:.1f}%)"


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1  Load base files
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 1: Load lifestyle_baseline_ptau.csv")

baseline = pd.read_csv(BASELINE_PTAU)
baseline["RID"] = baseline["RID"].astype(int)
N = len(baseline)
bl_rids = set(baseline["RID"])
print(f"  {N} patients, {baseline.shape[1]} cols")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2  Study entry diagnosis (dx_entry)
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 2: Study entry diagnosis — entry_research_group")

se = pd.read_csv(STUDY_ENTRY_FILE)
# Extract integer RID from subject_id format "002_S_0295"
se["RID"] = se["subject_id"].str.extract(r"_S_(\d+)$")[0].astype(int)
se = se.drop_duplicates("RID", keep="first")[["RID", "entry_research_group"]].rename(
    columns={"entry_research_group": "dx_entry"}
)

n_match = se[se["RID"].isin(bl_rids)]["RID"].nunique()
print(f"  Unique RIDs in study entry: {len(se)}")
print(f"  Matched to lifestyle_baseline: {pct(n_match, N)}")
print(f"  Value counts: {se['dx_entry'].value_counts().to_dict()}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3  DXSUM diagnosis (dx_bl)
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 3: DXSUM diagnosis — baseline (or latest visit fallback)")

dx = pd.read_csv(DXSUM_FILE)
dx["RID"] = pd.to_numeric(dx["RID"], errors="coerce")
dx["DIAGNOSIS"] = pd.to_numeric(dx["DIAGNOSIS"], errors="coerce")
dx = dx.dropna(subset=["RID", "DIAGNOSIS"])
dx["RID"] = dx["RID"].astype(int)

# Baseline rows
bl_dx = (
    dx[dx["VISCODE2"] == "bl"]
    .drop_duplicates("RID", keep="first")
    [["RID", "DIAGNOSIS", "VISCODE2"]]
)

# Fallback: for patients without baseline, use the latest available visit
missing_rids = bl_rids - set(bl_dx["RID"])
if missing_rids:
    fallback = (
        dx[dx["RID"].isin(missing_rids)]
        .assign(_date=lambda d: pd.to_datetime(d["EXAMDATE"], errors="coerce"))
        .sort_values(["RID", "_date"], ascending=[True, False])
        .drop_duplicates("RID", keep="first")
        [["RID", "DIAGNOSIS", "VISCODE2"]]
    )
    n_fb = len(fallback)
    print(
        f"  ⚠ {n_fb} patients had no baseline DXSUM → latest visit used as fallback")
    combined_dx = pd.concat([bl_dx, fallback], ignore_index=True)
else:
    combined_dx = bl_dx

combined_dx["dx_bl"] = combined_dx["DIAGNOSIS"].map(DX_MAP)
combined_dx["dx_bl_visit"] = combined_dx["VISCODE2"]
combined_dx = combined_dx[["RID", "dx_bl", "dx_bl_visit"]]

n_dx_match = combined_dx[combined_dx["RID"].isin(bl_rids)]["RID"].nunique()
n_bl_only = bl_dx[bl_dx["RID"].isin(bl_rids)]["RID"].nunique()
print(f"  Baseline DXSUM matched:        {pct(n_bl_only, N)}")
print(f"  Total (baseline + fallback):   {pct(n_dx_match, N)}")
dx_val_counts = (
    combined_dx[combined_dx["RID"].isin(bl_rids)]["dx_bl"]
    .value_counts().to_dict()
)
print(f"  Diagnosis breakdown: {dx_val_counts}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4  APOE genotype
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 4: APOE genotype — APOERES GENOTYPE")

ap = pd.read_csv(APOERES_FILE)
ap["RID"] = pd.to_numeric(ap["RID"], errors="coerce").astype("Int64")
apoe = (
    ap.dropna(subset=["RID", "GENOTYPE"])
    .drop_duplicates("RID", keep="first")
    [["RID", "GENOTYPE"]]
    .rename(columns={"GENOTYPE": "apoe_genotype"})
)
apoe["RID"] = apoe["RID"].astype(int)

n_apoe_match = apoe[apoe["RID"].isin(bl_rids)]["RID"].nunique()
print(f"  Unique RIDs in APOERES: {len(apoe)}")
print(f"  Matched to lifestyle_baseline: {pct(n_apoe_match, N)}")
print(
    f"  Genotype breakdown: {apoe['apoe_genotype'].value_counts().to_dict()}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5  Merge all → overwrite lifestyle_baseline_ptau.csv
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 5: Merge → lifestyle_baseline_ptau.csv")

# Drop existing columns if re-running
baseline = baseline.drop(
    columns=[c for c in NEW_COLS if c in baseline.columns])

result = (
    baseline
    .merge(se,           on="RID", how="left")
    .merge(combined_dx,  on="RID", how="left")
    .merge(apoe,         on="RID", how="left")
)

print(f"  Output: {len(result)} rows × {result.shape[1]} cols")
for col in NEW_COLS:
    n = result[col].notna().sum()
    print(f"  {col:18s}: {pct(n, N)}")

result.to_csv(BASELINE_PTAU, index=False)
print(f"  Saved: {BASELINE_PTAU}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6  Same merge for lifestyle_7var_ptau.csv
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 6: Merge → lifestyle_7var_ptau.csv")

var7 = pd.read_csv(VAR7_PTAU)
var7["RID"] = var7["RID"].astype(int)
var7 = var7.drop(columns=[c for c in NEW_COLS if c in var7.columns])

var7_out = (
    var7
    .merge(se,          on="RID", how="left")
    .merge(combined_dx, on="RID", how="left")
    .merge(apoe,        on="RID", how="left")
)

N7 = len(var7_out)
print(f"  Output: {N7} rows × {var7_out.shape[1]} cols")
for col in NEW_COLS:
    n = var7_out[col].notna().sum()
    print(f"  {col:18s}: {pct(n, N7)}")

var7_out.to_csv(VAR7_PTAU, index=False)
print(f"  Saved: {VAR7_PTAU}")

print("\n" + "═" * 72)
print("DONE — add_dx_apoe.py")
print("═" * 72)
