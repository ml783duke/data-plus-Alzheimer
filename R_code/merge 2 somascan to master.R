library(readr)
library(readxl)
library(dplyr)
library(writexl)

# File paths (relative to project root)
txt_file <- "protein_raw_data/somascan with gene ID.txt"

excel_file <- "protein_raw_data/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.xlsx"

output_file <- "protein_raw_data/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620_merged_tau.xlsx"

# Read txt file
txt_df <- read_csv(txt_file, show_col_types = FALSE)

# Read Excel file
protein_df <- read_excel(
  excel_file,
  sheet = "CruchagaLab_CSF_SOMAscan7k_Prot"
)

cat("Number of unique RIDs in txt SOMAscan file:", n_distinct(txt_df$RID), "\n")
cat("Number of unique RIDs in Excel/protein matrix SOMAscan file:", n_distinct(protein_df$RID), "\n")

# Variables to merge from txt file
vars_to_merge <- c(
  "RID",
  "ABETA40",
  "ABETA42",
  "TAU",
  "PTAU",
  "ptau181/ab42",
  "DX",
  "AGE",
  "PTGENDER",
  "PTEDUCAT",
  "PTETHCAT",
  "APOE4",
  "FDG",
  "CDRSB",
  "MMSE"
)

# Check whether all selected variables exist in txt file
missing_vars <- setdiff(vars_to_merge, names(txt_df))

if (length(missing_vars) > 0) {
  stop(
    "These variables are missing from the txt file: ",
    paste(missing_vars, collapse = ", ")
  )
}

# Keep one row per RID from txt file
subject_vars <- txt_df %>%
  select(all_of(vars_to_merge)) %>%
  group_by(RID) %>%
  summarise(
    across(everything(), ~ first(na.omit(.))),
    .groups = "drop"
  )

# Optional: check whether each RID has consistent values in txt file
consistency_check <- txt_df %>%
  select(all_of(vars_to_merge)) %>%
  group_by(RID) %>%
  summarise(
    across(
      everything(),
      ~ n_distinct(na.omit(.))
    ),
    .groups = "drop"
  ) %>%
  filter(if_any(-RID, ~ . > 1))

cat("Number of RIDs with inconsistent values:", nrow(consistency_check), "\n")

# Merge by RID
merged_df <- protein_df %>%
  left_join(subject_vars, by = "RID")

# Put newly merged variables after PlateId
new_vars <- setdiff(vars_to_merge, "RID")

merged_df <- merged_df %>%
  relocate(
    all_of(new_vars),
    .after = PlateId
  )

# Save as a new Excel file
write_xlsx(merged_df, output_file)

cat("Merged file saved to:\n", output_file, "\n")