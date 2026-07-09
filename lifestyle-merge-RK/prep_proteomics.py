"""
prep_proteomics.py — SOMAscan Proteomics Preprocessing
=======================================================

Inputs (protein_raw_data/):
  protein dict.csv                                       — analyte → gene/protein name lookup
  CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_...  — the actual protein expression matrix

Outputs (all in output/):
  SOMAscan7k_target.csv   — baseline, one row per patient, columns = protein Target names
  SOMAscan7k_gene.csv     — baseline, one row per patient, columns = EntrezGeneSymbol

Column naming:
  - Analyte maps to a unique Target/Gene → column name = that name (e.g. APOE, CLU).
  - Multiple aptamers share the same name → append SeqId: APOE.12345-6, APOE.67890-1.
    This keeps all aptamers in the file without silent column collisions.
  - 15 analytes have no EntrezGeneSymbol → fall back to Target name in gene file.

QC tips:
  - Columns with a "." (e.g. APOE.12345-6) are multi-aptamer proteins — check protein dict.csv
    if you need to know which aptamer is which.
  - Row count should be 708 baseline patients.
"""

import os
import pandas as pd

BASE = os.path.abspath(os.path.dirname(__file__))
PROTEIN_DIR = os.path.join(BASE, "protein_raw_data")
OUT_DIR = os.path.join(BASE, "output")
os.makedirs(OUT_DIR, exist_ok=True)

META_COLS = {"RID", "EXAMDATE", "GUSPECID", "Somalogic_Barcode_A",
             "VISCODE2", "ExtIdentifier", "PlateId"}

PROTEIN_MATRIX = os.path.join(
    PROTEIN_DIR,
    "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"
)
PROTEIN_DICT = os.path.join(PROTEIN_DIR, "protein dict.csv")


def section(title: str, note: str | None = None) -> None:
    print("\n" + "─" * 72)
    print(title)
    if note:
        print(note)
    print("─" * 72)


def make_unique_names(df: pd.DataFrame, name_col: str, seq_col: str = "SeqId") -> pd.Series:
    """
    Return a column-name Series from name_col, made unique by appending .SeqId
    when the same name appears for more than one aptamer.

    Example: if APOE appears for 3 aptamers → APOE.12345-6, APOE.67890-1, etc.
    """
    counts = df[name_col].value_counts()
    return df.apply(
        lambda r: r[name_col] if counts[r[name_col]] == 1
        else f"{r[name_col]}.{r[seq_col]}",
        axis=1
    )


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1  Load protein dict and build analyte → renamed-column maps
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 1: Load protein annotation dict (protein_raw_data/protein dict.csv)")

ann = pd.read_csv(PROTEIN_DICT)
proteins = ann[ann["Type"] == "Protein"].copy()
print(f"  Protein-type analytes: {len(proteins)}")
print(f"  Unique Target names:   {proteins['Target'].nunique()} "
      f"({(proteins['Target'].value_counts() > 1).sum()} with >1 aptamer)")
print(f"  Unique gene symbols:   {proteins['EntrezGeneSymbol'].nunique()} "
      f"({(proteins['EntrezGeneSymbol'].value_counts() > 1).sum()} with >1 aptamer)")
print(
    f"  Missing gene symbol:   {proteins['EntrezGeneSymbol'].isna().sum()} → fall back to Target")

# Fill missing gene symbols with Target name before building gene map
proteins["gene_fill"] = proteins["EntrezGeneSymbol"].fillna(proteins["Target"])

proteins["target_colname"] = make_unique_names(proteins, "Target")
proteins["gene_colname"] = make_unique_names(proteins, "gene_fill")

target_map = dict(zip(proteins["Analytes"], proteins["target_colname"]))
gene_map = dict(zip(proteins["Analytes"], proteins["gene_colname"]))

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2  Load protein matrix, filter to baseline, one row per patient
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 2: Load protein matrix → baseline, one row per patient",
    f"Source: {os.path.basename(PROTEIN_MATRIX)}",
)

prot_raw = pd.read_csv(PROTEIN_MATRIX)
print(f"  Full matrix:       {len(prot_raw)} rows × {prot_raw.shape[1]} cols")

prot_bl = (
    prot_raw[prot_raw["VISCODE2"] == "bl"]
    .copy()
    .drop_duplicates(subset="RID")
)
prot_bl["RID"] = prot_bl["RID"].astype(int)

analyte_cols = [c for c in prot_bl.columns if c not in META_COLS]
print(f"  Baseline patients: {len(prot_bl)}")
print(f"  Analyte columns:   {len(analyte_cols)}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3  Create SOMAscan7k_target.csv — columns named by protein Target
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 3: Create SOMAscan7k_target.csv")

target_analytes = [c for c in analyte_cols if c in target_map]
target_df = prot_bl[["RID", "VISCODE2"] +
                    target_analytes].rename(columns=target_map)

n_disambig_t = sum(1 for c in target_df.columns if "." in c)
print(f"  Analytes renamed: {len(target_analytes)} / {len(analyte_cols)}")
print(f"  Multi-aptamer cols (Name.SeqId): {n_disambig_t}")

target_out = os.path.join(OUT_DIR, "SOMAscan7k_target.csv")
target_df.to_csv(target_out, index=False)
print(f"  Saved: {target_out}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4  Create SOMAscan7k_gene.csv — columns named by EntrezGeneSymbol
# ═══════════════════════════════════════════════════════════════════════════════
section("STEP 4: Create SOMAscan7k_gene.csv")

gene_analytes = [c for c in analyte_cols if c in gene_map]
gene_df = prot_bl[["RID", "VISCODE2"] + gene_analytes].rename(columns=gene_map)

n_disambig_g = sum(1 for c in gene_df.columns if "." in c)
n_fallback = int(proteins["EntrezGeneSymbol"].isna().sum())
print(f"  Analytes renamed: {len(gene_analytes)} / {len(analyte_cols)}")
print(f"  Multi-aptamer cols (Gene.SeqId): {n_disambig_g}")
print(f"  Fallback to Target (no gene symbol): {n_fallback}")

gene_out = os.path.join(OUT_DIR, "SOMAscan7k_gene.csv")
gene_df.to_csv(gene_out, index=False)
print(f"  Saved: {gene_out}")

print("\n" + "═" * 72)
print("DONE — prep_proteomics.py")
print(
    f"  output/SOMAscan7k_target.csv  — {len(target_df)} patients × Target column names")
print(
    f"  output/SOMAscan7k_gene.csv    — {len(gene_df)} patients × gene symbol column names")
print("Next step: run merge.py")
print("═" * 72)
