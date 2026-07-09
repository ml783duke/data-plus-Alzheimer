###############################################################################
# fix_f06_forwardfill.R
# Fix f06_medhist variables in master_data:
# 1. Fix type conversion bug: read_excel auto-converts numeric columns to
#    logical when many values are 0/1, destroying actual numbers.
#    Solution: read as text, then convert to numeric (with -4 → NA).
# 2. Apply same-visit forward-fill: all 6 variables from the same visit.
#    Anchor: first visit where MH16SMOK is non-NA.
###############################################################################

library(readxl)
library(dplyr)

cat("========== Fixing f06 variables: correct types + forward-fill ==========\n")

# ---- 1. Read f06 with correct types ----
cat("\n--- Reading f06_medhist (as text to prevent type coercion) ---\n")

# Read ALL columns as text first to prevent auto type conversion
f06 <- read_excel("lifestyle_raw_data/f06_medhist.xlsx", col_types = "text")
cat(sprintf("Raw: %d rows x %d cols\n", nrow(f06), ncol(f06)))

med_cols <- c("MH14ALCH", "MH14AALCH", "MH16SMOK", "MH16ASMOK", "MH16BSMOK", "MH16CSMOK")

# Convert text to numeric: -4 is ADNI missing code
for (v in med_cols) {
  x <- suppressWarnings(as.numeric(f06[[v]]))
  x[x == -4] <- NA  # ADNI code for "not assessed"
  f06[[v]] <- x
}

cat("Converted to numeric (type coercion prevented, -4 → NA)\n")

# ---- 2. Visit ordering for f06 ----
# f06 visits: sc, f, m12, m24, m36, m48, m60, m72, m84
visit_order_f06 <- function(vcodes) {
  mapping <- c("sc" = -1, "f" = 0, "m12" = 12, "m24" = 24, "m36" = 36,
               "m48" = 48, "m60" = 60, "m72" = 72, "m84" = 84)
  order <- integer(length(vcodes))
  for (i in seq_along(vcodes)) {
    v <- vcodes[i]
    if (is.na(v)) {
      order[i] <- 9999
    } else if (v %in% names(mapping)) {
      order[i] <- mapping[v]
    } else {
      order[i] <- 9999
    }
  }
  return(order)
}

f06$visit_order <- visit_order_f06(f06$VISCODE2)
f06 <- f06[!is.na(f06$RID) & f06$visit_order < 9999, ]
f06 <- f06[order(f06$RID, f06$visit_order), ]

# ---- 3. Forward-fill with same-visit rule ----
# Anchor: first visit where MH16SMOK is non-NA
# Take ALL 6 variables from that same visit

all_rids <- unique(f06$RID)
cat(sprintf("Unique RIDs: %d\n", length(all_rids)))

result <- data.frame(RID = all_rids, stringsAsFactors = FALSE)
for (v in med_cols) result[[v]] <- NA_real_
result$MEDHIST_source_visit <- NA_character_

for (r in seq_along(all_rids)) {
  rows <- f06[f06$RID == all_rids[r], ]

  # Find first visit where the anchor (MH16SMOK) is non-NA
  anchor_idx <- which(!is.na(rows$MH16SMOK))[1]

  if (length(anchor_idx) == 1) {
    row <- rows[anchor_idx, ]
    for (v in med_cols) {
      result[r, v] <- row[[v]]
    }
    result$MEDHIST_source_visit[r] <- row$VISCODE2
  }
}

# Report
cat("\nSource visit distribution:\n")
print(table(result$MEDHIST_source_visit, useNA = "ifany"))

cat(sprintf("\nRIDs with data from non-sc visit: %d\n",
            sum(!is.na(result$MEDHIST_source_visit) & result$MEDHIST_source_visit != "sc")))

# Verify data types and values
cat("\nValue verification (numeric ranges):\n")
for (v in med_cols) {
  x <- na.omit(result[[v]])
  cat(sprintf("  %-12s: n=%d, range=[%s, %s], unique=%d\n",
              v, length(x),
              if(length(x)>0) min(x) else "NA",
              if(length(x)>0) max(x) else "NA",
              length(unique(x))))
}

# ---- 4. Load master and replace old f06 columns ----
cat("\n--- Updating master_data ---\n")
master <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Master before: %d rows x %d cols\n", nrow(master), ncol(master)))

# Remove old f06 columns
old_med_cols <- intersect(colnames(master), med_cols)
cat(sprintf("Removing %d old medhist columns\n", length(old_med_cols)))
master <- master[, setdiff(colnames(master), old_med_cols)]

# Merge new
master <- merge(master, result[, c("RID", med_cols)], by = "RID", all.x = TRUE)
cat(sprintf("After merge: %d rows x %d cols\n", nrow(master), ncol(master)))

# ---- 5. Reorder: medhist cols in lifestyle block before protein ----
is_protein <- grepl("^X[0-9]+\\.[0-9]+$", colnames(master))
prot_start <- which(is_protein)[1]
cat(sprintf("First protein at position: %d\n", prot_start))

cols_before <- setdiff(colnames(master)[1:(prot_start - 1)], med_cols)
cols_from_prot <- setdiff(colnames(master)[prot_start:length(colnames(master))], med_cols)
master <- master[, c(cols_before, med_cols, cols_from_prot)]

# Verify
new_prot_start <- which(grepl("^X[0-9]+\\.[0-9]+$", colnames(master)))[1]
med_pos <- which(colnames(master) %in% med_cols)
cat(sprintf("Medhist cols at %d-%d, protein at %d\n",
            min(med_pos), max(med_pos), new_prot_start))

# ---- 6. Report final coverage in master (baseline) ----
cat("\nMedhist coverage in master (baseline):\n")
df_bl <- master[master$VISCODE2 == "bl", ]
for (v in med_cols) {
  x <- na.omit(df_bl[[v]])
  cat(sprintf("  %-12s: n=%d/%d (%.1f%%), unique=%d, range=[%.2f, %.2f]\n",
              v, nrow(df_bl) - sum(is.na(df_bl[[v]])), nrow(df_bl),
              100 * (nrow(df_bl) - sum(is.na(df_bl[[v]]))) / nrow(df_bl),
              length(unique(x)),
              if(length(x)>0) min(x) else NA,
              if(length(x)>0) max(x) else NA))
}

# ---- 7. Save ----
writexl::write_xlsx(master, "master_data.xlsx")
cat(sprintf("\nSaved: master_data.xlsx (%d rows x %d cols)\n",
            nrow(master), ncol(master)))
cat("========== Done ==========\n")
