# =============================================================================
# ADNI Data Merge & Clean - Shared Reference & Utility Functions
# =============================================================================
# Source this file at the start of each analysis script:
#   source("R_code/00_reference.R")
# =============================================================================
# Last Updated: 2026-06-18

# ---- Libraries (install if missing) ----
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- Paths ---- (all short names — disk files already renamed)
DATA_DIR <- "."

FILE_SOMASCAN  <- file.path(DATA_DIR, "protein_raw_data", "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx")
FILE_DUKE_UPLC <- file.path(DATA_DIR, "lifestyle_raw_data", "f01_duke_uplc.xlsx")
FILE_LEIDEN_HPH <- file.path(DATA_DIR, "lifestyle_raw_data", "f02a_leiden_hph.xlsx")
FILE_LEIDEN_LPH <- file.path(DATA_DIR, "lifestyle_raw_data", "f02b_leiden_lph.xlsx")
FILE_DUKE_METAB <- file.path(DATA_DIR, "lifestyle_raw_data", "f03_duke_metabolon.xlsx")
FILE_NPI       <- file.path(DATA_DIR, "lifestyle_raw_data", "f04_npi.xlsx")
FILE_NPIQ      <- file.path(DATA_DIR, "lifestyle_raw_data", "f05_npiq.xlsx")
FILE_MEDHIST   <- file.path(DATA_DIR, "lifestyle_raw_data", "f06_medhist.xlsx")
FILE_BHR_BL    <- file.path(DATA_DIR, "lifestyle_raw_data", "f07_bhr_baseline.xlsx")
FILE_BHR_LONG  <- file.path(DATA_DIR, "lifestyle_raw_data", "f08_bhr_longitudinal.xlsx")
FILE_BHR_CG    <- file.path(DATA_DIR, "lifestyle_raw_data", "f09_bhr_caregiver.xlsx")
FILE_PTDEMOG   <- file.path(DATA_DIR, "lifestyle_raw_data", "f10_ptdemog.xlsx")
FILE_LIFESTYLE <- file.path(DATA_DIR, "lifestyle_raw_data", "lifestyle_dictionary.xlsx")

# ---- Short names map (from lifestyle_dictionary.xlsx) ----
FILE_SHORT_NAMES <- list(
  f01_duke_uplc        = FILE_DUKE_UPLC,
  f02a_leiden_hph      = FILE_LEIDEN_HPH,
  f02b_leiden_lph      = FILE_LEIDEN_LPH,
  f03_duke_metabolon   = FILE_DUKE_METAB,
  f04_npi              = FILE_NPI,
  f05_npiq             = FILE_NPIQ,
  f06_medhist          = FILE_MEDHIST,
  f07_bhr_baseline     = FILE_BHR_BL,
  f08_bhr_longitudinal = FILE_BHR_LONG,
  f09_bhr_caregiver    = FILE_BHR_CG,
  f10_ptdemog          = FILE_PTDEMOG
)

# ---- ID Conversion Functions ----

#' Convert PTID to RID
#'
#' Extracts the last 4 characters of PTID (the numeric ID after the final underscore).
#' Example: "027_S_4919" -> "4919", "027_S_0069" -> "0069"
#'
#' @param ptid Character vector of PTID identifiers
#' @return Character vector of RID identifiers (as character to preserve leading zeros)
ptid_to_rid <- function(ptid) {
  # Extract everything after the last underscore
  rid <- stringr::str_extract(ptid, "[^_]+$")
  # Ensure character type to preserve leading zeros
  as.character(rid)
}

#' Standardize RID for matching
#'
#' ✅ DECIDED: Pad to 4 digits with leading zeros.
#' Example: "69" -> "0069", "4919" -> "4919"
#'
#' @param rid Character or numeric vector of RIDs
#' @return Character vector of 4-digit zero-padded RIDs
rid_standardize <- function(rid) {
  # Pad to 4 digits with leading zeros
  stringr::str_pad(as.character(rid), width = 4, side = "left", pad = "0")
}

# ---- (Legacy alias — using the decided strategy) ----
rid_pad_4 <- rid_standardize

# ---- File Reading Helpers ----

#' Read an ADNI Excel file and standardize RID column
#'
#' Automatically detects which ID columns exist (RID, PTID, AboutRID, AboutPTID)
#' and ensures a clean, 4-digit-padded RID column is present.
#' If only PTID exists, converts to RID.
#'
#' @param filepath Path to Excel file
#' @return data.frame with a standardized 4-digit RID column
read_adni_excel <- function(filepath) {
  df <- readxl::read_excel(filepath)

  # Rename non-standard ID columns to standard names
  if ("AboutRID" %in% names(df)) {
    df <- df %>% rename(RID = AboutRID)
  }
  if ("AboutPTID" %in% names(df)) {
    df <- df %>% rename(PTID = AboutPTID)
  }

  # If RID missing but PTID exists, convert PTID to RID
  if (!"RID" %in% names(df) && "PTID" %in% names(df)) {
    df$RID <- ptid_to_rid(df$PTID)
    message("RID generated from PTID for: ", basename(filepath))
  }

  # Standardize RID to 4-digit character
  if ("RID" %in% names(df)) {
    df$RID <- rid_standardize(df$RID)
  }

  return(df)
}

#' Standardize time column to VISCODE2
#'
#' Handles files where time is named differently (e.g., BHR files use "Timepoint").
#' Renames to VISCODE2 and normalizes values to lowercase.
#'
#' **If no recognized time column is found, STOPS with a clear warning**
#' so the user can decide how to handle it manually.
#'
#' @param df data.frame
#' @param filename Optional file name for diagnostic messages
#' @return data.frame with standardized VISCODE2 column
standardize_time_col <- function(df, filename = "") {
  # Check for recognized time column names
  known_time_cols <- c("VISCODE2", "VISCODE", "Timepoint", "VISIT")

  found <- intersect(known_time_cols, names(df))

  if (length(found) == 0) {
    stop(
      "\n========================================\n",
      "WARNING: No recognized time column found in: ", filename, "\n",
      "Available columns: ", paste(head(names(df), 20), collapse = ", "), "...\n",
      "Expected one of: ", paste(known_time_cols, collapse = ", "), "\n",
      "Please check this file manually and decide how to proceed.\n",
      "========================================"
    )
  }

  # Prioritize VISCODE2, then VISCODE, then Timepoint
  if ("VISCODE2" %in% found) {
    # Already standard — just normalize
  } else if ("VISCODE" %in% found) {
    df <- df %>% rename(VISCODE2 = VISCODE)
    message("Renamed 'VISCODE' -> 'VISCODE2' for: ", filename)
  } else if ("Timepoint" %in% found) {
    df <- df %>% rename(VISCODE2 = Timepoint)
    message("Renamed 'Timepoint' -> 'VISCODE2' for: ", filename)
  }

  # Normalize VISCODE2 values to lowercase
  if ("VISCODE2" %in% names(df)) {
    df$VISCODE2 <- tolower(as.character(df$VISCODE2))
  }

  return(df)
}

#' Full join an ADNI table onto the SOMAscan anchor
#'
#' Always uses RID + VISCODE2 as the merge key to prevent incorrect
#' Cartesian expansion for longitudinal data.
#'
#' **Prerequisite**: Both anchor and new_df must already have
#' standardized RID (4-digit char) and VISCODE2 columns.
#' Use read_adni_excel() + standardize_time_col() first.
#'
#' @param anchor The SOMAscan main table (or accumulated master table)
#' @param new_df The new table to merge (must have RID + VISCODE2)
#' @param table_name Optional label for diagnostic messages
#' @return data.frame after full join on RID + VISCODE2
merge_adni <- function(anchor, new_df, table_name = "") {
  # Pre-check: both tables must have RID and VISCODE2
  missing_in_anchor <- setdiff(c("RID", "VISCODE2"), names(anchor))
  missing_in_new    <- setdiff(c("RID", "VISCODE2"), names(new_df))

  if (length(missing_in_anchor) > 0 || length(missing_in_new) > 0) {
    stop(
      "\n========================================\n",
      "WARNING: Cannot merge — missing required columns.\n",
      if (length(missing_in_anchor) > 0)
        paste0("  Anchor missing: ", paste(missing_in_anchor, collapse = ", "), "\n"),
      if (length(missing_in_new) > 0)
        paste0("  ", table_name, " missing: ", paste(missing_in_new, collapse = ", "), "\n"),
      "Make sure to run read_adni_excel() AND standardize_time_col() on both tables first.\n",
      "========================================"
    )
  }

  # Full join on RID + VISCODE2
  result <- dplyr::full_join(anchor, new_df, by = c("RID", "VISCODE2"))

  n_new_rids <- sum(!unique(new_df$RID) %in% unique(anchor$RID))
  n_new_rows <- nrow(result) - nrow(anchor)

  message(sprintf("Merged %s: +%d new RIDs, +%d total rows",
                  table_name, n_new_rids, n_new_rows))

  return(result)
}

# ---- Message ----
message("ADNI reference functions loaded. Data directory: ", normalizePath(DATA_DIR))
message("RID standardization: 4-digit padding (e.g., '69' -> '0069')")
message("Merge key: RID + VISCODE2 (dual-key match)")
message("Join strategy: FULL JOIN (no data loss from either side)")
