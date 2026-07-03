"""
prep_lifestyle.py — Lifestyle Data Preprocessing
=================================================

CONFIG (top of file):
  SKIP_WIDE_EXPANSION — set True to skip the slow Step 1 once it's been QC'd.

Outputs (all in output/):
  lifestyle_wide.csv                         — all visits, one row per patient, all sources merged
  deduplicated_raw/{source}_deduplicated.csv — per-source version of the above
  filtered/{source}_filtered.csv            — baseline-only, selected vars, per source
  lifestyle_baseline.csv                    — all sources merged, one row per patient at baseline

Column {src}_visit in lifestyle_baseline.csv:
  "bl"                  — row came from the target baseline visit
  "m06 (no bl)"         — patient had no baseline; earliest available visit used instead
  "earliest"            — file has no visit code (BHR files); earliest date row used
"""

import os
import re
import pandas as pd
from functools import reduce

BASE          = os.path.abspath(os.path.dirname(__file__))
LIFESTYLE_DIR = os.path.join(BASE, "lifestyle_raw_data")
OUT_DIR       = os.path.join(BASE, "output")
FILT_DIR      = os.path.join(OUT_DIR, "filtered")
DEDUP_DIR     = os.path.join(OUT_DIR, "deduplicated_raw")
os.makedirs(OUT_DIR,   exist_ok=True)
os.makedirs(FILT_DIR,  exist_ok=True)
os.makedirs(DEDUP_DIR, exist_ok=True)

# ── Configuration ─────────────────────────────────────────────────────────────
# Set True once you've QC'd deduplicated_raw/ — Step 1 is slow and only needs
# to be rerun if source files change.
SKIP_WIDE_EXPANSION = False
# ──────────────────────────────────────────────────────────────────────────────

FILE_RULES = {
    # filename: (viscode_col, viscode_val, rid_col, sort_col)
    # viscode_col/val — filter to this visit for baseline; None = no visit filter
    # rid_col         — patient ID column (renamed to RID on load)
    # sort_col        — date column used for chronological ordering
    "f01_duke_uplc.xlsx":        ("VISCODE2",  "bl",     "RID",            None),
    "f02a_leiden_hph.xlsx":      ("VISCODE2",  "bl",     "RID",            None),
    "f02b_leiden_lph.xlsx":      ("VISCODE2",  "bl",     "RID",            None),
    "f03_duke_metabolon.xlsx":   ("VISCODE2",  "bl",     "RID",            None),
    "f04_npi.xlsx":              ("VISCODE2",  "bl",     "RID",            None),
    "f05_npiq.xlsx":             ("VISCODE2",  "bl",     "RID",            None),
    "f06_medhist.xlsx":          ("VISCODE2",  "sc",     "RID",            None),
    "f07_bhr_baseline.xlsx":     (None,        None,     "RID",            "CollectedDate_DRVD"),
    "f08_bhr_longitudinal.xlsx": (None,        None,     "RID",            "CollectedDate_DRVD"),
    "f09_bhr_caregiver.xlsx":    ("Timepoint", "sp-m00", "StudyPartnerID", "CollectedDate_DRVD"),
    "f10_ptdemog.xlsx":          ("VISCODE2",  "sc",     "RID",            None),
}

DATE_COLS  = ["EXAMDATE", "CollectedDate_DRVD", "CollectedDateTime", "USERDATE2", "USERDATE"]
VISIT_COLS = ["VISCODE2", "VISCODE", "Timepoint"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_vars(var_str: str) -> list[str]:
    """Expand 'NPIK1-NPIK8; NPIK9A' into ['NPIK1', ..., 'NPIK8', 'NPIK9A']."""
    result = []
    for token in [v.strip() for v in str(var_str).split(";")]:
        m = re.match(r"^([A-Za-z_]+)(\d+)-[A-Za-z_]+(\d+)$", token)
        if m:
            prefix, n1, n2 = m.group(1), int(m.group(2)), int(m.group(3))
            result.extend(f"{prefix}{i}" for i in range(n1, n2 + 1))
        else:
            result.append(token)
    return result


def first_per_rid(df: pd.DataFrame, sort_col: str | None = None) -> pd.DataFrame:
    """Keep the first row per patient, optionally sorted by a date column first."""
    if sort_col and sort_col in df.columns:
        df = df.sort_values(sort_col)
    return df.drop_duplicates(subset="RID", keep="first")


def select_visit(
    df: pd.DataFrame,
    viscode_col: str | None,
    viscode_val: str | None,
    sort_col: str | None,
    flag_col: str,
) -> pd.DataFrame:
    """
    Return one row per patient at the target visit (viscode_val).

    For patients who have no row with that visit code, fall back to their
    earliest available row (sorted by date if possible) and mark the flag_col
    with the actual visit used + "(no {viscode_val})" so it's easy to spot.

    For files with no viscode filter (BHR files), take the earliest row per
    patient by date and mark flag_col as "earliest".
    """
    if not (viscode_col and viscode_col in df.columns and viscode_val):
        # No visit filter: sort by date and take earliest per patient
        date_col = sort_col if (sort_col and sort_col in df.columns) else next(
            (c for c in DATE_COLS if c in df.columns), None
        )
        result = first_per_rid(df, date_col)
        result = result.copy()
        result[flag_col] = "earliest"
        return result

    # Patients with the target visit
    target = df[df[viscode_col] == viscode_val].copy()
    target = first_per_rid(target, sort_col)  # deduplicate within target visit
    target[flag_col] = viscode_val

    # Patients with NO target visit → use earliest other row
    missing_rids = df[~df["RID"].isin(set(target["RID"]))]
    if len(missing_rids) == 0:
        return target

    # Sort fallback rows by date if possible
    date_col = sort_col if (sort_col and sort_col in missing_rids.columns) else next(
        (c for c in DATE_COLS if c in missing_rids.columns), None
    )
    if date_col:
        missing_rids = missing_rids.copy()
        missing_rids["_d"] = pd.to_datetime(missing_rids[date_col], errors="coerce")
        missing_rids = missing_rids.sort_values(["RID", "_d"]).drop(columns="_d")

    fallback = missing_rids.drop_duplicates("RID", keep="first").copy()
    fallback[flag_col] = fallback[viscode_col].astype(str) + f" (no {viscode_val})"
    print(f"    ⚠ {len(fallback)} patients had no '{viscode_val}' → earliest visit used instead")

    return pd.concat([target, fallback], ignore_index=True)


def expand_to_wide(df: pd.DataFrame, viscode_col: str | None, sort_col: str | None) -> pd.DataFrame:
    """
    Convert a source file that may have multiple rows per patient into one wide row.
    DHA with 2 visits → DHA_1, DHA_2.
    Adds metadata: n_raw_rows, extra_duplicate_rows, has_duplicate_rows,
                   visit_time_values, first/last observed date.
    """
    df = df.copy()
    df["_row_order"] = range(len(df))

    date_col  = sort_col if (sort_col and sort_col in df.columns) else next(
        (c for c in DATE_COLS if c in df.columns), None
    )
    visit_col = next((c for c in VISIT_COLS if c in df.columns), None)

    sort_keys = ["RID"]
    if date_col:
        df["_date"] = pd.to_datetime(df[date_col], errors="coerce")
        sort_keys  += ["_date", "_row_order"]
    elif viscode_col and viscode_col in df.columns:
        sort_keys  += [viscode_col, "_row_order"]
    else:
        sort_keys  += ["_row_order"]
    df = df.sort_values(sort_keys, kind="mergesort", na_position="last")

    df["_entry"] = df.groupby("RID").cumcount() + 1

    g    = df.groupby("RID", sort=False)
    meta = g.size().rename("n_raw_rows").reset_index()
    meta["extra_duplicate_rows"] = (meta["n_raw_rows"] - 1).clip(lower=0)
    meta["has_duplicate_rows"]   = meta["n_raw_rows"] > 1

    if visit_col:
        meta["visit_time_values"] = (
            g[visit_col]
            .apply(lambda s: "|".join(dict.fromkeys(s.dropna().astype(str))))
            .values
        )
    if date_col:
        date_range = g["_date"].agg(["min", "max"]).reset_index()
        date_range["first_observed_date"] = date_range["min"].dt.date.astype("string")
        date_range["last_observed_date"]  = date_range["max"].dt.date.astype("string")
        meta = meta.merge(
            date_range[["RID", "first_observed_date", "last_observed_date"]],
            on="RID", how="left"
        )

    skip       = {"RID", "_entry", "_row_order", "_date"}
    value_cols = [c for c in df.columns if c not in skip]
    wide       = df.pivot(index="RID", columns="_entry", values=value_cols)
    wide.columns = [f"{col}_{n}" for col, n in wide.columns]

    return meta.merge(wide.reset_index(), on="RID", how="left")


def prefix_source(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Rename non-RID columns: DHA_1 → DHA_f01_1, n_raw_rows → n_raw_rows_f01."""
    def rename(col: str) -> str:
        if col == "RID":
            return col
        m = re.match(r"^(.+)_(\d+)$", col)
        return f"{m.group(1)}_{prefix}_{m.group(2)}" if m else f"{col}_{prefix}"
    return df.rename(columns=rename)


def section(title: str, note: str | None = None) -> None:
    print("\n" + "─" * 72)
    print(title)
    if note:
        print(note)
    print("─" * 72)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1  Expand every source file → one wide row per patient (optional)
# ═══════════════════════════════════════════════════════════════════════════════
if SKIP_WIDE_EXPANSION:
    section("STEP 1: Wide expansion SKIPPED (SKIP_WIDE_EXPANSION = True)")
else:
    section(
        "STEP 1: Expand source files to wide format",
        "Each patient gets one row; repeated visits become DHA_1, DHA_2, etc.",
    )

    per_patient_sources: list[pd.DataFrame] = []

    for fname, (viscode_col, viscode_val, rid_col, sort_col) in FILE_RULES.items():
        fpath = os.path.join(LIFESTYLE_DIR, fname)
        if not os.path.exists(fpath):
            print(f"  SKIP (not found): {fname}")
            continue

        df = pd.read_excel(fpath)
        if rid_col != "RID":
            df = df.rename(columns={rid_col: "RID"})
        df["RID"] = pd.to_numeric(df["RID"], errors="coerce")
        df = df.dropna(subset=["RID"]).copy()
        df["RID"] = df["RID"].astype(int)

        wide       = expand_to_wide(df, viscode_col, sort_col)
        n_dup      = int((wide["n_raw_rows"] > 1).sum())
        n_extra    = int(wide["extra_duplicate_rows"].sum())
        src_prefix = fname.split("_")[0]

        dedup_path = os.path.join(DEDUP_DIR, f"{os.path.splitext(fname)[0]}_deduplicated.csv")
        wide.to_csv(dedup_path, index=False)

        per_patient_sources.append(prefix_source(wide, src_prefix))
        print(f"  {fname}: {len(wide)} patients | {n_dup} with >1 row | {n_extra} extra rows")
        print(f"    → {dedup_path}")

    section("Merging all sources on RID...")
    lifestyle_wide = reduce(
        lambda l, r: pd.merge(l, r, on="RID", how="outer"),
        per_patient_sources,
    ).sort_values("RID").reset_index(drop=True)

    print(f"  Shape: {lifestyle_wide.shape[0]} patients × {lifestyle_wide.shape[1]} cols")
    wide_out = os.path.join(OUT_DIR, "lifestyle_wide.csv")
    lifestyle_wide.to_csv(wide_out, index=False)
    print(f"  Saved: {wide_out}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2  Parse lifestyle_dictionary.xlsx → variable list per source file
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 2: Parse lifestyle dictionary",
    "Map each source file to the variables listed in lifestyle_dictionary.xlsx.",
)

dict_df = pd.read_excel(os.path.join(LIFESTYLE_DIR, "lifestyle_dictionary.xlsx"))
print(f"  Dictionary rows: {len(dict_df)}")

file_vars: dict[str, list[str]] = {}

for _, row in dict_df.iterrows():
    files = [f.strip() for f in str(row["New file name"]).split(";")]
    vars_ = parse_vars(row["Variable(s) in table"])

    if len(files) == 1:
        file_vars.setdefault(files[0], []).extend(vars_)
    else:
        for fname in files:
            fpath = os.path.join(LIFESTYLE_DIR, fname)
            if not os.path.exists(fpath):
                continue
            file_cols = pd.read_excel(fpath, nrows=0).columns.tolist()
            for v in vars_:
                if v in file_cols:
                    file_vars.setdefault(fname, []).append(v)

for fname in file_vars:
    seen: set[str] = set()
    file_vars[fname] = [v for v in file_vars[fname] if not (v in seen or seen.add(v))]  # type: ignore[func-returns-value]

print("\n  Variables per source file:")
for fname, vs in file_vars.items():
    print(f"    {fname}: {vs}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3  Load each source, select baseline (or earliest fallback) → lifestyle_baseline.csv
# ═══════════════════════════════════════════════════════════════════════════════
section(
    "STEP 3: Build lifestyle_baseline.csv",
    "One row per patient per source. Target visit used when available; "
    "earliest visit used as fallback (flagged in {src}_visit column).",
)

extracted: dict[str, pd.DataFrame] = {}

for fname, vars_ in file_vars.items():
    fpath = os.path.join(LIFESTYLE_DIR, fname)
    if not os.path.exists(fpath):
        print(f"  SKIP (not found): {fname}")
        continue

    viscode_col, viscode_val, rid_col, sort_col = FILE_RULES[fname]
    src_prefix = fname.split("_")[0]
    flag_col   = f"{src_prefix}_visit"

    # Load ALL rows (not pre-filtered); select_visit handles the baseline logic
    need     = [rid_col] + vars_
    if viscode_col: need = [viscode_col] + need
    if sort_col:    need = [sort_col] + need
    # Include any date columns present (needed for fallback chronological sorting)
    all_cols = pd.read_excel(fpath, nrows=0).columns.tolist()
    extra_dates = [c for c in DATE_COLS if c in all_cols and c not in need]
    df = pd.read_excel(fpath, usecols=[c for c in need + extra_dates if c in all_cols])

    if rid_col != "RID":
        df = df.rename(columns={rid_col: "RID"})
    df["RID"] = pd.to_numeric(df["RID"], errors="coerce")
    df = df.dropna(subset=["RID"]).copy()
    df["RID"] = df["RID"].astype(int)

    df = select_visit(df, viscode_col, viscode_val, sort_col, flag_col)

    for v in vars_:
        if v in df.columns:
            df[v] = pd.to_numeric(df[v], errors="coerce")

    # f05_npiq and f04_npi both contain NPIK; rename f05's copy to avoid collision
    if fname == "f05_npiq.xlsx":
        df = df.rename(columns={"NPIK": "NPIK_f05", "NPIKSEV": "NPIKSEV_f05"})

    present_vars = [v for v in vars_ if v in df.columns]
    df = df[["RID", flag_col] + present_vars]

    n_fallback = (df[flag_col] != viscode_val).sum() if viscode_val else 0
    print(f"  {fname}: {len(df)} patients | vars={present_vars}")
    if n_fallback:
        print(f"    ↳ {n_fallback} used earliest fallback — check '{flag_col}' column")

    extracted[fname] = df

    filt_path = os.path.join(FILT_DIR, f"{os.path.splitext(fname)[0]}_filtered.csv")
    df.to_csv(filt_path, index=False)
    print(f"    → {filt_path}")

section("Merging lifestyle sources on RID...")
lifestyle_baseline = reduce(
    lambda l, r: pd.merge(l, r, on="RID", how="outer"),
    extracted.values(),
).sort_values("RID").reset_index(drop=True)

print(f"  lifestyle_baseline: {lifestyle_baseline.shape[0]} rows × {lifestyle_baseline.shape[1]} cols")
baseline_out = os.path.join(OUT_DIR, "lifestyle_baseline.csv")
lifestyle_baseline.to_csv(baseline_out, index=False)
print(f"  Saved: {baseline_out}")

print("\n" + "═" * 72)
print("DONE — prep_lifestyle.py")
skip_note = " (Step 1 skipped)" if SKIP_WIDE_EXPANSION else ""
print(f"  output/lifestyle_baseline.csv{skip_note}")
if not SKIP_WIDE_EXPANSION:
    print("  output/lifestyle_wide.csv")
    print("  output/deduplicated_raw/  (11 files)")
print(f"  output/filtered/  (11 files)")
print("Next step: run merge.py")
print("═" * 72)
