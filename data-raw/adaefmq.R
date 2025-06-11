# Generate ADAEFMQ dataset

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(forcats)

# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADAEFMQ dataset
gen_adaefmq <- function() {
  ## create fmq DATA
  # Step 1: Read Excel sheets
  consolidated_list <- readRDS("source_data/FDA_FMQ_Consolidated_List.rds")
  fmq_references <- readRDS("source_data/FDA_FMQ_References.rds")

  # Step 2: Sort data
  consolidated_list01 <- consolidated_list %>%
    select(FMQ, PT, Final_Classification) %>%
    arrange(FMQ, PT)

  fmq_references01 <- fmq_references %>%
    rename(FMQ = FMQ_NAME, SOC = SYSTEM_ORGAN_CLASS) %>%
    arrange(FMQ)

  # Step 3: Merge datasets
  fmq01 <- left_join(consolidated_list01, fmq_references01, by = "FMQ")

  # Step 4: Sort and keep specific columns
  fmq02 <- fmq01 %>%
    select(PT, FMQ, SOC, Final_Classification) %>%
    arrange(PT, FMQ, SOC)

  # Create a counter for each fmq per pt
  fmq02_fmq_index <- fmq02 %>%
    select(-c("SOC", "Final_Classification")) %>%
    group_by(PT) %>%
    mutate(fmq_index = row_number()) %>%
    ungroup()

  # Step 2: Pivot the FMQ dataset to wide format
  transposed_fmq <- fmq02_fmq_index %>%
    pivot_wider(
      names_from = fmq_index,
      names_prefix = "fmq_",
      values_from = FMQ
    )

  # Create a counter for each fmq per soc
  fmq02_soc_index <- fmq02 %>%
    select(-c("FMQ", "Final_Classification")) %>%
    group_by(PT) %>%
    mutate(fmq_index = row_number()) %>%
    ungroup()

  # Step 2: Pivot the FMQ dataset to wide format
  transposed_soc <- fmq02_soc_index %>%
    pivot_wider(
      names_from = fmq_index,
      names_prefix = "soc_",
      values_from = SOC
    )

  # Create a counter for each fmq per soc
  fmq02_scope_index <- fmq02 %>%
    select(-c("FMQ", "SOC")) %>%
    group_by(PT) %>%
    mutate(fmq_index = row_number()) %>%
    ungroup()

  # Step 2: Pivot the FMQ dataset to wide format
  transposed_scope <- fmq02_scope_index %>%
    pivot_wider(
      names_from = fmq_index,
      names_prefix = "scope_",
      values_from = Final_Classification
    )

  # Step 6: Merge transposed datasets
  final_dataset <- transposed_fmq %>%
    left_join(transposed_soc, by = "PT") %>%
    left_join(transposed_scope, by = "PT")

  # Step 7: Process final dataset
  tfmq <- final_dataset %>%
    mutate(
      across(starts_with("fmq_"), as.character), # Ensure character type
      across(starts_with("soc_"), as.character),
      across(starts_with("scope_"), as.character),
      PT = toupper(PT)
    ) %>%
    rename(AEDECOD = PT)

  # Renaming variables and dropping unused ones
  names(tfmq) <- str_replace_all(names(tfmq), "fmq_", "FMQ") %>%
    str_replace_all("soc_", "FMQSOC") %>%
    str_replace_all("scope_", "SCOPE")

  ## now create adaefmq

  # source data
  source(file.path("data-raw", "adae.R"))

  gen <- adae %>%
    left_join(tfmq, by = "AEDECOD") %>% # Changed to left_join to keep all ADAE records
    rename_with(~ gsub("^FMQSOC", "SOCFMQ", .), starts_with("FMQSOC")) %>%
    # Pivot FMQ, SOCFMQ, and SCOPE columns simultaneously
    pivot_longer(
      cols = matches("^(FMQ|SOCFMQ|SCOPE)\\d+$"), # Select columns like FMQ1, SOCFMQ1, SCOPE1, etc.
      names_to = c(".value", "FMQ_NUM"), # Create columns based on the prefix (FMQ, SOCFMQ, SCOPE)
      names_pattern = "^(FMQ|SOCFMQ|SCOPE)(\\d+)$", # Define how to split the column names
      values_drop_na = FALSE # Keep rows even if some values (SOCFMQ, SCOPE) are NA initially
    ) %>%
    filter(!is.na(FMQ)) %>% # Keep only rows where the original FMQ value was not NA
    rename(
      FMQNAM = FMQ, # Rename the columns created by pivot_longer
      FMQSOC = SOCFMQ,
      FMQCLASS = SCOPE
    ) %>%
    select(-FMQ_NUM) # Remove the helper column containing the number suffix

  return(gen)
}

adaefmq <- gen_adaefmq()
