"""
Tie participant_subtype_assignments.csv to master_data.xlsx by participant ID,
using the lifestyle variable selection defined in MergedLifestyle.qmd.

Input:
  - /Users/reinakobayashi/Data+/data/participant_subtype_assignments.csv
      "Participant Unique ID" -- this is the ADNI RID (matches master_data.xlsx's RID column)
  - master_data.xlsx (in this directory)
      RID is the participant identifier; one row per RID at baseline (VISCODE2 == 'bl')

Lifestyle variables: the exact set selected in MergedLifestyle.qmd's
"select-respective-variables" + "merge-lifestyle" chunks (f01, f02a, f02b, f03,
f04, f06, f07, f09, f10 -- f05/NPIQ and f08/longitudinal are commented out
in that file's merge step, so excluded here too).

Output: prints match report; saves matched table with subtype + lifestyle
columns to results/subtype_lifestyle_matched.csv
"""

import os
import pandas as pd

BASE_DIR = "/Users/reinakobayashi/Data+/data-plus-Alzheimer"
SUBTYPE_FILE = "/Users/reinakobayashi/Data+/data/participant_subtype_assignments.csv"
MASTER_FILE = os.path.join(BASE_DIR, "master_data.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

# Lifestyle variables exactly as selected in MergedLifestyle.qmd's merge-lifestyle
# chunk (f05/NPIQ and f08/bhr_longitudinal are commented out there, so excluded here)
LIFESTYLE_VARS = [
    # f01_duke_uplc -- nutrition
    "DHA", "EPA", "HCys",
    # f02a_leiden_hph -- oxylipins, high confidence
    "BSH_FA20.5_w3", "BSH_FA22.6_w3",
    # f02b_leiden_lph -- oxylipins, low confidence
    "BSL_FA20.5_w3", "BSL_FA22.6_w3",
    # f03_duke_metabolon -- fatty acids / cotinine biomarkers
    "X424", "X439", "X100008930", "X2050", "X100000665",
    "X100001181", "X848", "X100002717", "X100002719", "X100004494",
    # f04_npi -- neuropsychiatric inventory (sleep items)
    "NPIK", "NPIKTOT", "NPIK1", "NPIK2", "NPIK3", "NPIK4", "NPIK5",
    "NPIK6", "NPIK7", "NPIK8", "NPIK9A", "NPIK9B", "NPIK9C",
    # f06_medhist -- smoking + alcohol
    "MH14ALCH", "MH14AALCH", "MH16SMOK", "MH16ASMOK", "MH16BSMOK", "MH16CSMOK",
    # f07_bhr_baseline -- social (reduced activities)
    "QID46_9",
    # f09_bhr_caregiver -- social (loneliness), joined on StudyPartnerID == RID
    "QID44_1_3", "QID44_2_3",
    # f10_ptdemog -- physical activity / occupation
    "PTWORK", "PTWORKHS",
]

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: Load subtype assignments
# ═══════════════════════════════════════════════════════════════════════════
print("Loading participant_subtype_assignments.csv...")
subtype = pd.read_csv(SUBTYPE_FILE).rename(columns={"Participant Unique ID": "RID"})
n_subtype = len(subtype)
print(f"  {n_subtype} participants with subtype assignments")
print(f"  RID range: {subtype['RID'].min()}-{subtype['RID'].max()}\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: Load master_data.xlsx, keep baseline rows only (RID is otherwise
# duplicated across visits), and pull RID + DX + lifestyle columns
# ═══════════════════════════════════════════════════════════════════════════
print("Loading master_data.xlsx (RID, VISCODE2, DX, lifestyle columns)...")
usecols = ["RID", "VISCODE2", "DX"] + LIFESTYLE_VARS
master = pd.read_excel(MASTER_FILE, sheet_name="Sheet1", usecols=usecols)
print(f"  Master: {len(master)} rows, {master['RID'].nunique()} unique RIDs")

master_bl = master[master["VISCODE2"] == "bl"].drop_duplicates(subset="RID")
print(f"  Master baseline rows: {len(master_bl)} (one per RID)\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: Merge subtype assignments to master on RID
# ═══════════════════════════════════════════════════════════════════════════
merged = subtype.merge(master_bl, on="RID", how="left", indicator=True)
n_matched = (merged["_merge"] == "both").sum()
n_unmatched = (merged["_merge"] == "left_only").sum()

print("=" * 60)
print("MATCH REPORT")
print("=" * 60)
print(f"Subtype assignment participants:  {n_subtype}")
print(f"Matched to master_data.xlsx (baseline): {n_matched} ({n_matched/n_subtype*100:.1f}%)")
print(f"Unmatched: {n_unmatched} ({n_unmatched/n_subtype*100:.1f}%)")

if n_unmatched:
    unmatched_rids = merged.loc[merged["_merge"] == "left_only", "RID"].tolist()
    print(f"\nUnmatched RIDs: {unmatched_rids}")

merged = merged.drop(columns="_merge")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: Lifestyle variable coverage among matched participants
# ═══════════════════════════════════════════════════════════════════════════
matched = merged[merged["RID"].isin(master_bl["RID"])]
print(f"\nLifestyle variable coverage among {len(matched)} matched participants:")
for v in LIFESTYLE_VARS:
    n_present = matched[v].notna().sum()
    print(f"  {v}: {n_present}/{len(matched)} ({n_present/len(matched)*100:.1f}%)")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: Save
# ═══════════════════════════════════════════════════════════════════════════
out_path = os.path.join(OUT_DIR, "subtype_lifestyle_matched.csv")
merged.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
