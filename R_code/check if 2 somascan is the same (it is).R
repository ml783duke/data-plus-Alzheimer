library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# File paths (relative to project root)
txt_file <- "protein_raw_data/somascan with gene ID.txt"
wide_file <- "protein_raw_data/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx"

# Read files
long_df <- read_csv(txt_file, show_col_types = FALSE)

wide_df <- read_excel(
  wide_file,
  sheet = "CruchagaLab_CSF_SOMAscan7k_Prot"
)

# Identify protein columns in the wide file
protein_cols_wide <- names(wide_df)[str_detect(names(wide_df), "^X[0-9]+\\.[0-9]+$")]

# Identify protein IDs in the long file
protein_ids_long <- sort(unique(long_df$`Protein ID`))
protein_ids_wide <- sort(protein_cols_wide)

# Check whether both files contain the same proteins
proteins_only_in_long <- setdiff(protein_ids_long, protein_ids_wide)
proteins_only_in_wide <- setdiff(protein_ids_wide, protein_ids_long)

cat("Number of proteins in long file:", length(protein_ids_long), "\n")
cat("Number of proteins in wide file:", length(protein_ids_wide), "\n")
cat("Proteins only in long file:", length(proteins_only_in_long), "\n")
cat("Proteins only in wide file:", length(proteins_only_in_wide), "\n")

# Check whether both files contain the same RIDs
rids_long <- sort(unique(long_df$RID))
rids_wide <- sort(unique(wide_df$RID))

rids_only_in_long <- setdiff(rids_long, rids_wide)
rids_only_in_wide <- setdiff(rids_wide, rids_long)

cat("Number of RIDs in long file:", length(rids_long), "\n")
cat("Number of RIDs in wide file:", length(rids_wide), "\n")
cat("RIDs only in long file:", length(rids_only_in_long), "\n")
cat("RIDs only in wide file:", length(rids_only_in_wide), "\n")

# Convert wide protein matrix to long format
wide_long_df <- wide_df %>%
  select(RID, all_of(protein_cols_wide)) %>%
  mutate(
    across(
      all_of(protein_cols_wide),
      ~ suppressWarnings(as.numeric(.))
    )
  ) %>%
  pivot_longer(
    cols = all_of(protein_cols_wide),
    names_to = "Protein ID",
    values_to = "wide_value"
  ) %>%
  mutate(
    wide_log2_value = log2(wide_value)
  )

# Keep protein measurements from the long file
long_protein_df <- long_df %>%
  select(RID, `Protein ID`, Y) %>%
  mutate(
    Y = suppressWarnings(as.numeric(Y))
  )

# Join the two versions
compare_df <- long_protein_df %>%
  inner_join(
    wide_long_df,
    by = c("RID", "Protein ID")
  ) %>%
  mutate(
    difference = Y - wide_log2_value,
    abs_difference = abs(difference),
    is_match = case_when(
      is.na(Y) & is.na(wide_log2_value) ~ TRUE,
      is.na(Y) | is.na(wide_log2_value) ~ FALSE,
      abs_difference < 1e-6 ~ TRUE,
      TRUE ~ FALSE
    )
  )

# Summary of measurement agreement
comparison_summary <- compare_df %>%
  summarise(
    n_compared = n(),
    n_match = sum(is_match, na.rm = TRUE),
    n_mismatch = sum(!is_match, na.rm = TRUE),
    max_abs_difference = max(abs_difference, na.rm = TRUE),
    mean_abs_difference = mean(abs_difference, na.rm = TRUE)
  )

print(comparison_summary)

# Show examples of mismatches, if any
mismatch_examples <- compare_df %>%
  filter(!is_match) %>%
  arrange(desc(abs_difference)) %>%
  select(RID, `Protein ID`, Y, wide_value, wide_log2_value, difference, abs_difference) %>%
  head(20)

print(mismatch_examples)

# Optional: check duplicate measurements in the long file
duplicate_long_measurements <- long_df %>%
  count(RID, `Protein ID`) %>%
  filter(n > 1)

cat("Duplicate RID-protein pairs in long file:", nrow(duplicate_long_measurements), "\n")

# Optional: check duplicate RIDs in the wide file
duplicate_wide_rids <- wide_df %>%
  count(RID) %>%
  filter(n > 1)

cat("Duplicate RIDs in wide file:", nrow(duplicate_wide_rids), "\n")