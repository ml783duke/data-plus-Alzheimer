###############################################################################
# merge_lifestyle_to_master.R
# Extract baseline lifestyle variables from source files and merge into master
# - Each source file: filter to baseline (VISCODE2 == "bl")
# - Extract collaborator-specified lifestyle columns
# - Join by RID → insert after MMSE, before protein columns
###############################################################################

# ===========================================================================
# 0. Setup
# ===========================================================================
library(readxl)
library(dplyr)

DATA_DIR <- "lifestyle_raw_data"

cat("========== Merging Lifestyle Variables into Master ==========\n")

# ===========================================================================
# 1. Read and extract baseline lifestyle from each source file
# ===========================================================================

# --- f01: Duke UPLC (Nutrition: DHA, EPA, Homocysteine) ---
cat("\n--- f01_duke_uplc ---\n")
f01 <- read_excel(file.path(DATA_DIR, "f01_duke_uplc.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f01), ncol(f01)))
f01_bl <- f01[f01$VISCODE2 == "bl", c("RID", "DHA", "EPA", "HCys")]
f01_bl <- f01_bl[!is.na(f01_bl$RID), ]
# Convert to numeric
for (v in c("DHA", "EPA", "HCys")) {
  f01_bl[[v]] <- suppressWarnings(as.numeric(f01_bl[[v]]))
}
cat(sprintf("  Baseline: %d rows, %d unique RIDs\n", nrow(f01_bl), length(unique(f01_bl$RID))))

# --- f02a: Leiden HPH (Oxylipins, high confidence) ---
cat("\n--- f02a_leiden_hph ---\n")
f02a <- read_excel(file.path(DATA_DIR, "f02a_leiden_hph.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f02a), ncol(f02a)))
f02a_bl <- f02a[f02a$VISCODE2 == "bl", c("RID", "BSH_FA20.5_w3", "BSH_FA22.6_w3")]
f02a_bl <- f02a_bl[!is.na(f02a_bl$RID), ]
f02a_bl$BSH_FA20.5_w3 <- suppressWarnings(as.numeric(f02a_bl$BSH_FA20.5_w3))
f02a_bl$BSH_FA22.6_w3 <- suppressWarnings(as.numeric(f02a_bl$BSH_FA22.6_w3))
cat(sprintf("  Baseline: %d rows, %d unique RIDs\n", nrow(f02a_bl), length(unique(f02a_bl$RID))))

# --- f02b: Leiden LPH (Oxylipins, low confidence) ---
cat("\n--- f02b_leiden_lph ---\n")
f02b <- read_excel(file.path(DATA_DIR, "f02b_leiden_lph.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f02b), ncol(f02b)))
f02b_bl <- f02b[f02b$VISCODE2 == "bl", c("RID", "BSL_FA20.5_w3", "BSL_FA22.6_w3")]
f02b_bl <- f02b_bl[!is.na(f02b_bl$RID), ]
f02b_bl$BSL_FA20.5_w3 <- suppressWarnings(as.numeric(f02b_bl$BSL_FA20.5_w3))
f02b_bl$BSL_FA22.6_w3 <- suppressWarnings(as.numeric(f02b_bl$BSL_FA22.6_w3))
cat(sprintf("  Baseline: %d rows, %d unique RIDs\n", nrow(f02b_bl), length(unique(f02b_bl$RID))))

# --- f03: Duke Metabolon (Fatty acids + cotinine biomarkers) ---
cat("\n--- f03_duke_metabolon ---\n")
f03 <- read_excel(file.path(DATA_DIR, "f03_duke_metabolon.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f03), ncol(f03)))
metab_cols <- c("RID", "X424", "X439", "X100008930", "X2050", "X100000665",
                "X100001181", "X848", "X100002717", "X100002719", "X100004494")
f03_bl <- f03[f03$VISCODE2 == "bl", metab_cols]
f03_bl <- f03_bl[!is.na(f03_bl$RID), ]
for (v in setdiff(metab_cols, "RID")) {
  f03_bl[[v]] <- suppressWarnings(as.numeric(f03_bl[[v]]))
}
cat(sprintf("  Baseline: %d rows, %d unique RIDs\n", nrow(f03_bl), length(unique(f03_bl$RID))))

# --- f04: NPI (Sleep items + total) ---
cat("\n--- f04_npi ---\n")
f04 <- read_excel(file.path(DATA_DIR, "f04_npi.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f04), ncol(f04)))
npi_cols <- c("RID", "NPIK", "NPIKTOT", "NPIK1", "NPIK2", "NPIK3", "NPIK4",
              "NPIK5", "NPIK6", "NPIK7", "NPIK8", "NPIK9A", "NPIK9B", "NPIK9C")
f04_bl <- f04[f04$VISCODE2 == "bl", npi_cols]
f04_bl <- f04_bl[!is.na(f04_bl$RID), ]
for (v in setdiff(npi_cols, "RID")) {
  f04_bl[[v]] <- suppressWarnings(as.numeric(f04_bl[[v]]))
}
cat(sprintf("  Baseline: %d rows, %d unique RIDs\n", nrow(f04_bl), length(unique(f04_bl$RID))))

# --- f05: NPIQ (Sleep items, severity) ---
cat("\n--- f05_npiq ---\n")
f05 <- read_excel(file.path(DATA_DIR, "f05_npiq.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f05), ncol(f05)))
f05_bl <- f05[f05$VISCODE2 == "bl", c("RID", "NPIK", "NPIKSEV")]
f05_bl <- f05_bl[!is.na(f05_bl$RID), ]
# Rename to avoid conflict with f04 NPIK
colnames(f05_bl)[2:3] <- c("NPIQ_K", "NPIQ_KSEV")
f05_bl$NPIQ_K <- suppressWarnings(as.numeric(f05_bl$NPIQ_K))
f05_bl$NPIQ_KSEV <- suppressWarnings(as.numeric(f05_bl$NPIQ_KSEV))
cat(sprintf("  Baseline: %d rows, %d unique RIDs\n", nrow(f05_bl), length(unique(f05_bl$RID))))

# --- f06: MEDHIST (Smoking + Alcohol) ---
# NOTE: f06 has no "bl" visit; uses "sc" (screening) as first visit
cat("\n--- f06_medhist ---\n")
f06 <- read_excel(file.path(DATA_DIR, "f06_medhist.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f06), ncol(f06)))
cat(sprintf("  VISCODE2 values: %s\n", paste(names(table(f06$VISCODE2)), collapse=", ")))
medhist_cols <- c("RID", "MH14ALCH", "MH14AALCH", "MH16SMOK",
                  "MH16ASMOK", "MH16BSMOK", "MH16CSMOK")
# Use "sc" (screening) — no baseline "bl" in this file
f06_bl <- f06[f06$VISCODE2 == "sc", medhist_cols]
f06_bl <- f06_bl[!is.na(f06_bl$RID), ]
for (v in setdiff(medhist_cols, "RID")) {
  f06_bl[[v]] <- suppressWarnings(as.numeric(f06_bl[[v]]))
}
# Keep first non-NA row per RID if duplicates exist
f06_bl <- f06_bl[!duplicated(f06_bl$RID), ]
cat(sprintf("  Screening (sc): %d rows, %d unique RIDs\n", nrow(f06_bl), length(unique(f06_bl$RID))))

# --- f07: BHR Baseline (Social: reduced activities) ---
cat("\n--- f07_bhr_baseline ---\n")
f07 <- read_excel(file.path(DATA_DIR, "f07_bhr_baseline.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f07), ncol(f07)))
# This file is already baseline only; no VISCODE2. Use Timepoint to confirm
if ("Timepoint" %in% colnames(f07)) {
  cat(sprintf("  Timepoint values: %s\n", paste(unique(f07$Timepoint), collapse=", ")))
}
f07_bl <- f07[, c("RID", "QID46_9")]
f07_bl <- f07_bl[!is.na(f07_bl$RID), ]
f07_bl$QID46_9 <- suppressWarnings(as.numeric(f07_bl$QID46_9))
cat(sprintf("  Kept: %d rows, %d unique RIDs\n", nrow(f07_bl), length(unique(f07_bl$RID))))

# --- f09: BHR Caregiver (Social: loneliness) ---
cat("\n--- f09_bhr_caregiver ---\n")
f09 <- read_excel(file.path(DATA_DIR, "f09_bhr_caregiver.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f09), ncol(f09)))
# Uses StudyPartnerID (same as RID). Single measurement.
f09_bl <- f09[, c("StudyPartnerID", "QID44_1_3", "QID44_2_3")]
colnames(f09_bl)[1] <- "RID"
f09_bl <- f09_bl[!is.na(f09_bl$RID), ]
f09_bl$QID44_1_3 <- suppressWarnings(as.numeric(f09_bl$QID44_1_3))
f09_bl$QID44_2_3 <- suppressWarnings(as.numeric(f09_bl$QID44_2_3))
cat(sprintf("  Kept: %d rows, %d unique RIDs\n", nrow(f09_bl), length(unique(f09_bl$RID))))

# --- f10: PTDEMOG (Physical activity: occupation) ---
# NOTE: f10 has very few "bl"; most data at "sc" (screening)
cat("\n--- f10_ptdemog ---\n")
f10 <- read_excel(file.path(DATA_DIR, "f10_ptdemog.xlsx"))
cat(sprintf("  Raw: %d rows x %d cols\n", nrow(f10), ncol(f10)))
cat(sprintf("  VISCODE2 values: %s\n", paste(names(table(f10$VISCODE2)), collapse=", ")))
f10_bl <- f10[f10$VISCODE2 == "sc", c("RID", "PTWORK", "PTWORKHS")]
f10_bl <- f10_bl[!is.na(f10_bl$RID), ]
f10_bl$PTWORK <- suppressWarnings(as.numeric(f10_bl$PTWORK))
f10_bl$PTWORKHS <- suppressWarnings(as.numeric(f10_bl$PTWORKHS))
# Keep first non-NA row per RID if duplicates exist
f10_bl <- f10_bl[!duplicated(f10_bl$RID), ]
cat(sprintf("  Screening (sc): %d rows, %d unique RIDs\n", nrow(f10_bl), length(unique(f10_bl$RID))))

# ===========================================================================
# 2. Merge all lifestyle data frames by RID (one row per RID)
# ===========================================================================
cat("\n========== Merging lifestyle sources ==========\n")

lifestyle_list <- list(f01_bl, f02a_bl, f02b_bl, f03_bl, f04_bl, f05_bl,
                       f06_bl, f07_bl, f09_bl, f10_bl)

# Deduplicate each source: keep first row per RID
for (i in seq_along(lifestyle_list)) {
  n_before <- nrow(lifestyle_list[[i]])
  lifestyle_list[[i]] <- lifestyle_list[[i]][!duplicated(lifestyle_list[[i]]$RID), ]
  n_after <- nrow(lifestyle_list[[i]])
  if (n_before != n_after) {
    cat(sprintf("  Dedup source %d: %d → %d rows\n", i, n_before, n_after))
  }
}

# Merge sequentially by RID
lifestyle_merged <- Reduce(function(x, y) {
  merge(x, y, by = "RID", all = TRUE)
}, lifestyle_list)

cat(sprintf("Merged lifestyle: %d rows x %d columns\n",
            nrow(lifestyle_merged), ncol(lifestyle_merged)))

# Report coverage
n_total <- nrow(lifestyle_merged)
for (v in setdiff(colnames(lifestyle_merged), "RID")) {
  n_present <- sum(!is.na(lifestyle_merged[[v]]))
  cat(sprintf("  %s: %d / %d (%.1f%%)\n", v, n_present, n_total,
              100 * n_present / n_total))
}

# ===========================================================================
# 3. Load master and merge lifestyle
# ===========================================================================
cat("\n========== Merging into master_data.xlsx ==========\n")

master <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Master: %d rows x %d columns\n", nrow(master), ncol(master)))

# Find MMSE position
mmse_pos <- which(colnames(master) == "MMSE")
cat(sprintf("MMSE is column %d\n", mmse_pos))

# Merge lifestyle by RID
master <- merge(master, lifestyle_merged, by = "RID", all.x = TRUE)
cat(sprintf("After merge: %d rows x %d columns\n", nrow(master), ncol(master)))

# ===========================================================================
# 4. Reorder: lifestyle columns after MMSE, before protein columns
# ===========================================================================
cat("\n========== Reordering columns ==========\n")

lifestyle_vars <- setdiff(colnames(lifestyle_merged), "RID")
all_cols <- colnames(master)

# Columns before and including MMSE
cols_before_mmse <- all_cols[1:mmse_pos]

# Protein columns (start with X followed by digits and dot)
protein_cols <- grep("^X[0-9]+\\.[0-9]+$", all_cols, value = TRUE)

# Other columns (between MMSE and proteins, plus after proteins)
cols_other <- setdiff(all_cols, c(cols_before_mmse, lifestyle_vars, protein_cols))

# Rebuild: [before MMSE] + [lifestyle] + [protein] + [other]
new_order <- c(cols_before_mmse, lifestyle_vars, protein_cols, cols_other)
master <- master[, new_order]

cat(sprintf("Lifestyle columns (%d) placed after MMSE (col %d)\n",
            length(lifestyle_vars), mmse_pos))
cat(sprintf("First lifestyle col: %s (now col %d)\n",
            lifestyle_vars[1], mmse_pos + 1))
cat(sprintf("First protein col: %s (now col %d)\n",
            protein_cols[1], mmse_pos + length(lifestyle_vars) + 1))

# ===========================================================================
# 5. Save
# ===========================================================================
cat("\n========== Saving ==========\n")
writexl::write_xlsx(master, "master_data.xlsx")
cat(sprintf("Saved master_data.xlsx: %d rows x %d columns\n",
            nrow(master), ncol(master)))

# ===========================================================================
# 6. Summary report
# ===========================================================================
cat("\n========== Summary ==========\n")
cat(sprintf("Master rows: %d\n", nrow(master)))
cat(sprintf("Lifestyle variables added: %d\n", length(lifestyle_vars)))
cat(sprintf("  Nutrition (DHA, EPA, HCys, oxylipins, fatty acids): %d vars\n",
            sum(c("DHA","EPA","HCys","BSH_FA20.5_w3","BSH_FA22.6_w3",
                  "BSL_FA20.5_w3","BSL_FA22.6_w3","X424","X439","X100008930",
                  "X2050","X100000665","X100001181") %in% lifestyle_vars)))
cat(sprintf("  Smoking (status, pack-years, cotinine): %d vars\n",
            sum(c("MH16SMOK","MH16ASMOK","MH16BSMOK","MH16CSMOK",
                  "X848","X100002717","X100002719","X100004494") %in% lifestyle_vars)))
cat(sprintf("  Alcohol (abuse, drinks/day): %d vars\n",
            sum(c("MH14ALCH","MH14AALCH") %in% lifestyle_vars)))
cat(sprintf("  Sleep (NPI + NPIQ): %d vars\n",
            sum(c("NPIK","NPIKTOT","NPIK1","NPIK2","NPIK3","NPIK4","NPIK5",
                  "NPIK6","NPIK7","NPIK8","NPIK9A","NPIK9B","NPIK9C",
                  "NPIQ_K","NPIQ_KSEV") %in% lifestyle_vars)))
cat(sprintf("  Social (activities, loneliness): %d vars\n",
            sum(c("QID46_9","QID44_1_3","QID44_2_3") %in% lifestyle_vars)))
cat(sprintf("  Physical activity (occupation): %d vars\n",
            sum(c("PTWORK","PTWORKHS") %in% lifestyle_vars)))

# RID match stats
master_rids <- master$RID[!is.na(master$RID)]
lifestyle_rids <- lifestyle_merged$RID
n_matched <- sum(master_rids %in% lifestyle_rids)
cat(sprintf("\nRID matching: %d / %d master RIDs have lifestyle data\n",
            n_matched, length(master_rids)))

cat("\n========== Done ==========\n")
