"""
merge.py — Merge Lifestyle + Proteomics
========================================

Requires (run these first):
  python3 prep_lifestyle.py    → output/master_lifestyle.csv
  python3 prep_proteomics.py   → output/SOMAscan7k_gene.csv

Outputs (all in output/):
  lifestyle_proteomics.csv   — inner join of lifestyle + SOMAscan7k_gene
                                      (protein columns named by EntrezGeneSymbol)
  lifestyle_7var.csv         — 7 key lifestyle variables only
"""

import os
import pandas as pd

BASE = os.path.abspath(os.path.dirname(__file__))
OUT_DIR = os.path.join(BASE, "output")

KEY_VARS = ["DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK"]


def section(title: str, note: str | None = None) -> None:
    print("\n" + "─" * 72)
    print(title)
    if note:
        print(note)
    print("─" * 72)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1  Load upstream outputs
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 1: Load lifestyle and proteomics tables")

lifestyle_path = os.path.join(OUT_DIR, "lifestyle_baseline.csv")
prot_path = os.path.join(OUT_DIR, "SOMAscan7k_gene.csv")

for path in [lifestyle_path, prot_path]:
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Missing: {path}\n"
            "Run prep_lifestyle.py and prep_proteomics.py first."
        )

master_lifestyle = pd.read_csv(lifestyle_path)
prot_gene = pd.read_csv(prot_path)

print(
    f"  master_lifestyle: {master_lifestyle.shape[0]} rows × {master_lifestyle.shape[1]} cols")
print(
    f"  SOMAscan7k_gene:  {prot_gene.shape[0]} rows × {prot_gene.shape[1]} cols")

n_overlap = master_lifestyle["RID"].isin(prot_gene["RID"]).sum()
print(f"  RID overlap: {n_overlap} lifestyle patients present in proteomics")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2  Inner join → lifestyle_proteomics.csv
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 2: Merge → lifestyle_proteomics.csv",
    "Inner join keeps only patients present in both datasets.",
)

protein_cols = [c for c in prot_gene.columns if c not in {"RID", "VISCODE2"}]
master_full = master_lifestyle.merge(
    prot_gene[["RID"] + protein_cols], on="RID", how="inner"
)
print(f"  Result: {len(master_full)} patients × {master_full.shape[1]} cols")
print(
    f"  Lifestyle-only patients dropped: {master_lifestyle['RID'].nunique() - len(master_full)}")
print(
    f"  Proteomics-only patients dropped: {prot_gene['RID'].nunique() - len(master_full)}")

full_out = os.path.join(OUT_DIR, "lifestyle_proteomics.csv")
master_full.to_csv(full_out, index=False)
print(f"  Saved: {full_out}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3  7-variable lifestyle subset → lifestyle_7var.csv
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 3: Create master_lifestyle_7var.csv",
    f"Subset to: {KEY_VARS}",
)

# Exclude visit flag columns ({src}_visit) — they're in lifestyle_baseline.csv but not KEY_VARS
present = [v for v in KEY_VARS if v in master_lifestyle.columns]
missing = [v for v in KEY_VARS if v not in master_lifestyle.columns]
if missing:
    print(f"  WARNING: not found in master_lifestyle: {missing}")

master_7var = (
    master_lifestyle[["RID"] + present]
    .dropna(subset=present, how="all")
)
print(f"  Patients with ≥1 key variable: {len(master_7var)}")
print("  Coverage:")
for v in present:
    n = master_7var[v].notna().sum()
    print(f"    {v}: {n}/{len(master_7var)} ({n/len(master_7var)*100:.1f}%)")

var7_out = os.path.join(OUT_DIR, "lifestyle_7var.csv")
master_7var.to_csv(var7_out, index=False)
print(f"  Saved: {var7_out}")

print("\n" + "═" * 72)
print("DONE — merge.py")
print(
    f"  output/lifestyle_proteomics.csv  — {len(master_full)} patients, gene-named protein cols")
print(
    f"  output/lifestyle_7var.csv         — {len(master_7var)} patients")
print("═" * 72)
