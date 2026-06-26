###############################################################################
# add_new_ptau_compare.R
# Add new_ptau column (from NULISA pTau-217 CSF NPQ) to master_data.xlsx
# - Match by RID
# - new_ptau inserted after PlateId column
# - Report matching statistics
###############################################################################

# ===========================================================================
# 0. Setup
# ===========================================================================
library(readxl)

# Paths (relative to project root)
master_path <- "master_data.xlsx"
nulisa_path <- "protein_raw_data/BSHRI_PLA_CSF_NULISA_CNS_22Jun2026.xlsx"
out_path    <- "master_data.xlsx"  # overwrite

# ===========================================================================
# 1. Load master data
# ===========================================================================
cat("\n========== Loading master_data ==========\n")
master <- read_excel(master_path, sheet = 1)
cat(sprintf("Master: %d rows x %d columns\n", nrow(master), ncol(master)))

# Find PlateId position
plate_pos <- which(colnames(master) == "PlateId")
cat(sprintf("PlateId is column %d\n", plate_pos))

# ===========================================================================
# 2. Load and filter NULISA data
# ===========================================================================
cat("\n========== Loading NULISA data ==========\n")
nulisa <- read_excel(nulisa_path, sheet = 1, col_types = "text")
cat(sprintf("NULISA raw: %d rows x %d columns\n", nrow(nulisa), ncol(nulisa)))

# Filter: Target == "pTau-217" & SampleMatrixType == "CSF" & RID not NA
nulisa_sub <- nulisa[
  nulisa$Target == "pTau-217" &
  nulisa$SampleMatrixType == "CSF" &
  !is.na(nulisa$RID),
]
cat(sprintf("After filtering (pTau-217 + CSF + has RID): %d rows\n", nrow(nulisa_sub)))

# Convert RID and NPQ to numeric
nulisa_sub$RID_num <- as.numeric(nulisa_sub$RID)
nulisa_sub$NPQ_num  <- as.numeric(nulisa_sub$NPQ)

# Remove rows where conversion failed
nulisa_sub <- nulisa_sub[!is.na(nulisa_sub$RID_num) & !is.na(nulisa_sub$NPQ_num), ]
cat(sprintf("After numeric conversion: %d rows, %d unique RIDs\n",
            nrow(nulisa_sub), length(unique(nulisa_sub$RID_num))))

# ===========================================================================
# 3. Handle duplicate RIDs (keep first occurrence)
# ===========================================================================
cat("\n========== Handling duplicate RIDs ==========\n")
dup_rids <- names(table(nulisa_sub$RID_num)[table(nulisa_sub$RID_num) > 1])
if (length(dup_rids) > 0) {
  cat(sprintf("Duplicate RIDs found: %d\n", length(dup_rids)))
  for (r in dup_rids) {
    rows <- which(nulisa_sub$RID_num == as.numeric(r))
    cat(sprintf("  RID %s: %d rows, NPQ = %s -> keeping first (%.4f)\n",
                r, length(rows),
                paste(round(nulisa_sub$NPQ_num[rows], 2), collapse = ", "),
                nulisa_sub$NPQ_num[rows[1]]))
  }
  # Keep first
  nulisa_sub <- nulisa_sub[!duplicated(nulisa_sub$RID_num), ]
  cat(sprintf("After dedup: %d rows\n", nrow(nulisa_sub)))
} else {
  cat("No duplicate RIDs found.\n")
}

# ===========================================================================
# 4. Match and add new_ptau
# ===========================================================================
cat("\n========== Matching RIDs ==========\n")

master_rids <- unique(master$RID[!is.na(master$RID)])
nulisa_rids <- unique(nulisa_sub$RID_num)

cat(sprintf("Master unique RIDs: %d\n", length(master_rids)))
cat(sprintf("NULISA pTau-217 CSF unique RIDs: %d\n", length(nulisa_rids)))

# Create lookup: RID -> NPQ
nulisa_lookup <- nulisa_sub[, c("RID_num", "NPQ_num")]
colnames(nulisa_lookup) <- c("RID", "new_ptau")

# Match
master$new_ptau <- nulisa_lookup$new_ptau[match(master$RID, nulisa_lookup$RID)]

n_matched <- sum(!is.na(master$new_ptau))
cat(sprintf("Matched (new_ptau added): %d / %d master rows\n", n_matched, nrow(master)))

# ===========================================================================
# 5. Insert new_ptau after PlateId
# ===========================================================================
cat("\n========== Reordering columns ==========\n")

# Remove new_ptau from end, re-insert after PlateId
new_ptau_col <- master$new_ptau
master$new_ptau <- NULL  # remove from end

# Rebuild column order
cols_before <- colnames(master)[1:plate_pos]
cols_after  <- colnames(master)[(plate_pos + 1):ncol(master)]

master$new_ptau <- new_ptau_col  # add at end
master <- master[, c(cols_before, "new_ptau", cols_after)]

cat(sprintf("new_ptau inserted at column %d (after PlateId)\n", plate_pos + 1))

# ===========================================================================
# 6. Report matches and mismatches
# ===========================================================================
cat("\n========== Matching Report ==========\n")

in_both      <- intersect(master_rids, nulisa_rids)
in_master_only <- setdiff(master_rids, nulisa_rids)
in_nulisa_only <- setdiff(nulisa_rids, master_rids)

cat(sprintf("RIDs in both:            %d\n", length(in_both)))
cat(sprintf("RIDs in master only:     %d (no pTau-217 CSF data)\n", length(in_master_only)))
cat(sprintf("RIDs in NULISA only:     %d (not in master, skipped)\n", length(in_nulisa_only)))

# Among master rows with PTAU, how many have new_ptau?
has_ptau <- master[!is.na(master$PTAU), ]
cat(sprintf("\nMaster rows with PTAU value: %d\n", nrow(has_ptau)))
cat(sprintf("  Of these, have new_ptau: %d\n", sum(!is.na(has_ptau$new_ptau))))

# Summary stats
cat("\n========== new_ptau Summary ==========\n")
cat(sprintf("new_ptau range: %.2f - %.2f\n", min(master$new_ptau, na.rm = TRUE),
            max(master$new_ptau, na.rm = TRUE)))
cat(sprintf("new_ptau missing: %d / %d (%.1f%%)\n",
            sum(is.na(master$new_ptau)), nrow(master),
            100 * sum(is.na(master$new_ptau)) / nrow(master)))

# Compare PTAU vs new_ptau
both_ptau <- master[!is.na(master$PTAU) & !is.na(master$new_ptau), ]
if (nrow(both_ptau) > 0) {
  cat(sprintf("\nRows with both PTAU and new_ptau: %d\n", nrow(both_ptau)))
  cat(sprintf("PTAU range:     %.2f - %.2f\n", min(both_ptau$PTAU), max(both_ptau$PTAU)))
  cat(sprintf("new_ptau range: %.2f - %.2f\n", min(both_ptau$new_ptau), max(both_ptau$new_ptau)))
  cor_val <- cor(both_ptau$PTAU, both_ptau$new_ptau, method = "spearman", use = "complete.obs")
  cat(sprintf("Spearman correlation (PTAU vs new_ptau): rho = %.4f\n", cor_val))
}

# ===========================================================================
# 7. Save
# ===========================================================================
cat("\n========== Saving ==========\n")
writexl::write_xlsx(master, out_path)
cat(sprintf("Saved: %s\n", out_path))
cat(sprintf("Final dimensions: %d rows x %d columns\n", nrow(master), ncol(master)))

cat("\n========== Done ==========\n")
