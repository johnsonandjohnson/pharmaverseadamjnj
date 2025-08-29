# Generate ADAGOCMQ dataset

# Load necessary libraries
library(dplyr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))


#' Derive Records Combining `ATERMN` Levels
#' 
#' For example, if one subject has one record with ATERMN = 121 and another
#' record with ATERMN = 122, and the start date of two those records is within
#' 7 days, then set ATERMN = 12, create a record. If there are multiple
#' combinations satisfy the above criteria, create one record for each
#' combination.
#' 
#' @param df A data frame.
#' @param levels A numeric vector of two `ATERMN` levels.
#' @param level A single numeric that will be assigned to the derived records.
#' 
#' @return The input data frame with the derived records appended.
#' 
#' @noRd
#' 
#' @examples
#' df <- dplyr::tibble(
#'   USUBJID = c('a', 'a', 'a'),
#'   ASTDT = as.Date(c("2020-11-01", "2020-11-02", "2020-11-03")),
#'   OCMQNAM = rep("Hypersensitivity", 3),
#'   ATERMN = c(121, 122, 121)
#' )
#' 
#' derive_combined_atermn(df, c(121, 122), 12)
derive_combined_atermn <- function(df, levels, level) {
  ids <- unique(df[["USUBJID"]])
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(ATERMN %in% levels)
    
    if (!all(levels %in% subject[["ATERMN"]])) next
    
    for (i in seq_len(NROW(subject))) {
      current <- subject[i, ]
      
      records_within_7 <- subject[-(1:i), ] |>
        dplyr::filter(abs(ASTDT - current[["ASTDT"]]) <= 7) |>
        dplyr::filter(ATERMN != current[["ATERMN"]]) |>
        dplyr::mutate(ATERMN = level) |>
        
        # Clear some variables to prevent derived records being processed again
        dplyr::mutate(dplyr::across(
          dplyr::any_of(c("OCMQNAM", "HYPSCAT")),
          \(x) NA_character_
        ))
      
      df <- dplyr::bind_rows(df, records_within_7)
    }
  }
  
  df
}


# Generate ADAGOCMQ dataset
gen_adagocmq <- function() {
  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    HYPSCAT = "Hypersensitivity Category"
  )
  
  # Create records related to ADAEOCMQ ----------------------------------------

  # Get source data
  raw_adaeocmq <- adaeocmq
  
  # For `HYPSCAT` derivation. Downloaded from:
  # https://www.fda.gov/drugs/development-resources/
  # office-new-drugs-custom-medical-queries-ocmqs
  terms <-
    file.path("data-raw", "OCMQs_v3.0.xlsm") |>
    readxl::read_excel(sheet = "Hypersensitivity") |>
    dplyr::mutate(
      AEDECOD = toupper(Term),
      HYPSCAT = substr(`Algorithmic Category`, 1, 1)
    ) |>
    dplyr::select(AEDECOD, HYPSCAT)
  
  gen_adaeocmq <- raw_adaeocmq |>
    # Derive `HYPSCAT`
    dplyr::left_join(terms) |>
    
    # Derive `ATERMN`
    dplyr::mutate(ATERMN = dplyr::case_when(
      OCMQNAM == "Hypersensitivity" & HYPSCAT == "A" ~ 11,
      OCMQNAM == "Hypersensitivity" & HYPSCAT == "B" ~ 121,
      OCMQNAM == "Hypersensitivity" & HYPSCAT == "C" ~ 122
    )) |>
    derive_combined_atermn(c(121, 122), 12) |>
    dplyr::mutate(ATERMN = dplyr::case_when(
      OCMQNAM == "Hypersensitivity" & HYPSCAT == "D" ~ 131,
      .default = ATERMN
    )) |>
    derive_combined_atermn(c(121, 131), 13) |>
    derive_combined_atermn(c(122, 131), 14) |>
    dplyr::mutate(ATERMN = dplyr::case_when(
      OCMQNAM == "Hyperglycemia" & OCMQCLSS == "Narrow" ~ 21,
      .default = ATERMN
    )) |>
    
    # Derive `ATERM`
    dplyr::mutate(ATERM = dplyr::case_when(
      ATERMN == 11 ~ "Any hypersensitivity OCMQ narrow term",
      ATERMN == 121 ~ "Respiratory",
      ATERMN == 122 ~ "Skin Reaction",
      ATERMN == 12 ~ "Respiratory + Skin Reaction",
      ATERMN == 131 ~ "Systemic Reaction",
      ATERMN == 13 ~ "Respiratory + Systemic Reaction",
      ATERMN == 14 ~ "Skin + Systemic Reaction",
      ATERMN == 21 ~ "Any Hyperglycemia OCMQ Narrow Term"
    )) |>
    
    # Handle NA values and convert characters to factors
    df_na(char_as_factor = TRUE) |>
    # Restore labels
    restore_labels(raw_adaeocmq, additional_labels)
  
  return(gen)
}


# Generate the dataset
adagocmq <- gen_adagocmq()
