###############################################################################
# fix_npi_forwardfill.R
# Fix NPI variables in master_data: forward-fill from baseline
# Key rule: ALL NPI variables for a given subject come from the SAME visit.
# For each RID, find the first visit where NPIK (anchor variable) is non-NA,
# then take ALL NPI variables from that exact visit.
# If NPIK is never available, fall back to the first visit with any NPI data.
###############################################################################

library(readxl)
library(dplyr)

cat("========== Fixing NPI variables: forward-fill (same-visit rule) ==========\n")

# ---- 1. Build forward-filled NPI data from f04 ----
cat("\n--- Reading f04_npi ---\n")
f04 <- read_excel("lifestyle_raw_data/f04_npi.xlsx")

build_visit_order <- function(vcodes) {
  order <- integer(length(vcodes))
  for (i in seq_along(vcodes)) {
    v <- vcodes[i]
    if (is.na(v)) {
      order[i] <- 9999
    } else if (v == "sc") {
      order[i] <- -1
    } else if (v == "bl") {
      order[i] <- 0
    } else {
      order[i] <- as.integer(gsub("m", "", v))
    }
  }
  return(order)
}

npi_cols <- c("NPIK", "NPIKTOT",
              paste0("NPIK", 1:8),
              "NPIK9A", "NPIK9B", "NPIK9C")

f04$visit_order <- build_visit_order(f04$VISCODE2)
f04_npi <- f04[f04$visit_order >= 0 & !is.na(f04$RID), ]
f04_npi <- f04_npi[order(f04_npi$RID, f04_npi$visit_order), ]

all_rids <- unique(f04_npi$RID)
cat(sprintf("Unique RIDs (bl+): %d\n", length(all_rids)))

result <- data.frame(RID = all_rids, stringsAsFactors = FALSE)
for (v in npi_cols) result[[v]] <- NA_real_
result$NPI_source_visit <- NA_character_

for (r in seq_along(all_rids)) {
  rows <- f04_npi[f04_npi$RID == all_rids[r], ]

  # Find the first visit where NPIK (anchor) is non-NA
  anchor_idx <- which(!is.na(rows$NPIK))[1]

  if (length(anchor_idx) == 1) {
    # Take ALL NPI variables from this visit
    row <- rows[anchor_idx, ]
    for (v in npi_cols) {
      if (!is.na(row[[v]])) result[r, v] <- as.numeric(row[[v]])
    }
    result$NPI_source_visit[r] <- row$VISCODE2
  } else {
    # NPIK never available: try NPIKTOT as fallback
    fallback_idx <- which(!is.na(rows$NPIKTOT))[1]
    if (length(fallback_idx) == 1) {
      row <- rows[fallback_idx, ]
      for (v in npi_cols) {
        if (!is.na(row[[v]])) result[r, v] <- as.numeric(row[[v]])
      }
      result$NPI_source_visit[r] <- row$VISCODE2
    }
  }
}

# Report visit distribution
cat(sprintf("Forward-filled NPI: %d RIDs, %d cols\n", nrow(result), ncol(result) - 1))
cat("\nSource visit distribution:\n")
print(table(result$NPI_source_visit, useNA = "ifany"))

cat(sprintf("\nRIDs with NPIK from non-baseline visit: %d\n",
            sum(!is.na(result$NPI_source_visit) & result$NPI_source_visit != "bl")))

# ---- 2. Load master and remove old NPI columns ----
cat("\n--- Updating master_data ---\n")
master <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Master before: %d rows x %d cols\n", nrow(master), ncol(master)))

old_npi_cols <- intersect(colnames(master), npi_cols)
cat(sprintf("Removing %d old NPI columns\n", length(old_npi_cols)))

master <- master[, setdiff(colnames(master), old_npi_cols)]

# ---- 3. Merge new NPI data ----
master <- merge(master, result, by = "RID", all.x = TRUE)
cat(sprintf("After merge: %d rows\n", nrow(master)))

# ---- 4. Reorder: NPI cols right before protein columns ----
# Identify protein columns
is_protein <- grepl("^X[0-9]+\\.[0-9]+$", colnames(master))
prot_start <- which(is_protein)[1]
cat(sprintf("First protein column at position %d: %s\n",
            prot_start, colnames(master)[prot_start]))

# Split columns
cols_before <- colnames(master)[1:(prot_start - 1)]
cols_before <- setdiff(cols_before, npi_cols)
cols_from_prot <- colnames(master)[prot_start:length(colnames(master))]
cols_from_prot <- setdiff(cols_from_prot, npi_cols)

new_order <- c(cols_before, npi_cols, cols_from_prot)
master <- master[, new_order]

# Verify
new_prot_start <- which(grepl("^X[0-9]+\\.[0-9]+$", colnames(master)))[1]
npi_pos <- which(colnames(master) %in% npi_cols)
cat(sprintf("NPI cols at %d-%d, protein starts at %d\n",
            min(npi_pos), max(npi_pos), new_prot_start))

# ---- 5. Report coverage in master (baseline) ----
cat("\nNPI coverage in master (baseline subjects):\n")
df_bl <- master[master$VISCODE2 == "bl", ]
for (v in npi_cols) {
  n <- sum(!is.na(df_bl[[v]]))
  cat(sprintf("  %-10s: %d / %d (%.1f%%)\n", v, n, nrow(df_bl), 100*n/nrow(df_bl)))
}

# ---- 6. Save ----
writexl::write_xlsx(master, "master_data.xlsx")
cat(sprintf("\nSaved: master_data.xlsx (%d rows x %d cols)\n",
            nrow(master), ncol(master)))
cat("========== Done ==========\n")
