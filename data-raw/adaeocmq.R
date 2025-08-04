# Generate ADAEOCMQ dataset

# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(forcats)

# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADAEOCMQ dataset
gen_adaeocmq <- function() {
  ## create ocmq DATA
  # Step 1: Read RAW Data
  consolidated_list <- readRDS("source_data/FDA_OCMQ_Consolidated_List.rds")
  ocmq_references <- readRDS("source_data/FDA_OCMQ_References.rds")

  # Step 2: Sort data
  consolidated_list01 <- consolidated_list %>%
    select(OCMQ, Term, Final_Classification = `Final Classification`) %>%
    mutate(Term = stringr::str_to_sentence(Term)) %>%
    arrange(OCMQ, Term)

  ocmq_references01 <- ocmq_references %>%
    rename(OCMQ = `OCMQ NAME`, SOC = `ORGAN SYSTEM`) %>%
    mutate(SOC = stringr::str_to_sentence(SOC)) %>%
    arrange(OCMQ)

  # Step 3: Merge datasets
  ocmq01 <- left_join(consolidated_list01, ocmq_references01, by = "OCMQ")

  # Step 4: Sort and keep specific columns
  ocmq02 <- ocmq01 %>%
    select(Term, OCMQ, SOC, Final_Classification) %>%
    arrange(Term, OCMQ, SOC)

  # Create a counter for each ocmq per pt
  ocmq02_ocmq_index <- ocmq02 %>%
    select(-c("SOC", "Final_Classification")) %>%
    group_by(Term) %>%
    mutate(ocmq_index = row_number()) %>%
    ungroup()

  # Step 2: Pivot the OCMQ dataset to wide format
  transposed_ocmq <- ocmq02_ocmq_index %>%
    pivot_wider(
      names_from = ocmq_index,
      names_prefix = "ocmq_",
      values_from = OCMQ
    )

  # Create a counter for each ocmq per soc
  ocmq02_soc_index <- ocmq02 %>%
    select(-c("OCMQ", "Final_Classification")) %>%
    group_by(Term) %>%
    mutate(ocmq_index = row_number()) %>%
    ungroup()

  # Step 2: Pivot the OCMQ dataset to wide format
  transposed_soc <- ocmq02_soc_index %>%
    pivot_wider(
      names_from = ocmq_index,
      names_prefix = "soc_",
      values_from = SOC
    )

  # Create a counter for each ocmq per soc
  ocmq02_scope_index <- ocmq02 %>%
    select(-c("OCMQ", "SOC")) %>%
    group_by(Term) %>%
    mutate(ocmq_index = row_number()) %>%
    ungroup()

  # Step 2: Pivot the OCMQ dataset to wide format
  transposed_scope <- ocmq02_scope_index %>%
    pivot_wider(
      names_from = ocmq_index,
      names_prefix = "scope_",
      values_from = Final_Classification
    )

  # Step 6: Merge transposed datasets
  final_dataset <- transposed_ocmq %>%
    left_join(transposed_soc, by = "Term") %>%
    left_join(transposed_scope, by = "Term")

  # Step 7: Process final dataset
  tocmq <- final_dataset %>%
    mutate(
      across(starts_with("ocmq_"), as.character), # Ensure character type
      across(starts_with("soc_"), as.character),
      across(starts_with("scope_"), as.character),
      Term = toupper(Term)
    ) %>%
    rename(AEDECOD = Term)

  # Renaming variables and dropping unused ones
  names(tocmq) <- str_replace_all(names(tocmq), "ocmq_", "OCMQ") %>%
    str_replace_all("soc_", "OCMQSOC") %>%
    str_replace_all("scope_", "SCOPE")

  ## now create adaeocmq

  # source data
  source(file.path("data-raw", "adae.R"))

  gen <- adae %>%
    left_join(tocmq, by = "AEDECOD") %>% # Changed to left_join to keep all ADAE records
    rename_with(~ gsub("^OCMQSOC", "SOCOCMQ", .), starts_with("OCMQSOC")) %>%
    # Pivot OCMQ, SOCOCMQ, and SCOPE columns simultaneously
    pivot_longer(
      cols = matches("^(OCMQ|SOCOCMQ|SCOPE)\\d+$"), # Select columns like OCMQ1, SOCOCMQ1, SCOPE1, etc.
      names_to = c(".value", "OCMQ_NUM"), # Create columns based on the prefix (OCMQ, SOCOCMQ, SCOPE)
      names_pattern = "^(OCMQ|SOCOCMQ|SCOPE)(\\d+)$", # Define how to split the column names
      values_drop_na = FALSE # Keep rows even if some values (SOCOCMQ, SCOPE) are NA initially
    ) %>%
    filter(!is.na(OCMQ)) %>%
    rename(
      OCMQNAM = OCMQ, # Rename the columns created by pivot_longer
      OCMQSOC = SOCOCMQ,
      OCMQCLSS = SCOPE
    ) %>%
    select(-OCMQ_NUM)

  # Add gender-specific flags
  gen <- gen %>%
    mutate(
      GENSPMFL = ifelse(OCMQNAM %in% c("Erectile Dysfunction", "Gynecomastia"), "Y", NA_character_),
      GENSPFFL = ifelse(OCMQNAM %in% c(
        "Abnormal Uterine Bleeding", "Amenorrhea", "Bacterial Vaginosis",
        "Decreased Menstrual Bleeding", "Excessive Menstrual Bleeding"
      ), "Y", NA_character_)
    )

  # Additional labels for new variables not in the source dataset
  additional_labels <- list(
    OCMQNAM  = "Custom Medical Query Name",
    OCMQSOC  = "Custom Medical Query System Organ Class",
    OCMQCLSS = "Custom Medical Query Scope",
    GENSPMFL = "Gender Specific OCMQ Male Flag",
    GENSPFFL = "Gender Specific OCMQ Female Flag"
  )

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels
  gen <- restore_labels(
    df = gen,
    orig_df = gen,
    additional_labels = additional_labels
  )

  return(gen)
}

adaeocmq <- gen_adaeocmq()
