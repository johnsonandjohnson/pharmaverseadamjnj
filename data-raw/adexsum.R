# Generate ADEXSUM dataset

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(formatters)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADEXSUM dataset
gen_adexsum <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Start with pharmaverseadam::adex for structure
  raw <- pharmaverseadam::adex

  # Get unique subject IDs from the base data
  unique_subjects <- unique(raw$USUBJID)

  # Define the specific parameter values needed
  param_mapping <- data.frame(
    PARAMCD = c(
      "CUMDOSE",
      "DOSEDAYS",
      "FINDD",
      "MEANDD",
      "MEANDDI",
      "MODEDD",
      "NUMCYC",
      "TNUMDOS",
      "TRTCOMP",
      "TRTDURM",
      "TRTDURY"
    ),
    PARAM = c(
      "Cumulative dose ([unit])",
      "Total dosing days of treatment (excluding days off treatment)",
      "Final daily dose ([unit]/day)",
      "Average daily dose ([unit]/day) (excluding days off treatment)",
      "Average daily dose ([unit]/day) (including days off treatment)",
      "Modal daily dose ([unit]/day)",
      "Total number of cycles received",
      "Total number of administrations",
      "Compliance (%)",
      "Duration of treatment, months",
      "Duration of treatment, years"
    ),
    stringsAsFactors = FALSE
  )

  # Create empty list to store records
  records <- list()
  record_index <- 1

  # For each subject, create one record for each parameter
  for (subj in unique_subjects) {
    for (i in seq_len(nrow(param_mapping))) {
      # Use the base structure but replace key fields
      record <- raw[1, ] # Use first row as template
      record$USUBJID <- subj
      record$PARAMCD <- param_mapping$PARAMCD[i]
      record$PARAM <- param_mapping$PARAM[i]

      # Generate appropriate AVAL based on parameter type
      record$AVAL <- dplyr::case_when(
        record$PARAMCD == "TRTDURM" ~ as.numeric(round(runif(1, 0, 38), 1)),
        record$PARAMCD == "TRTDURY" ~ as.numeric(round(runif(1, 0, 3.5), 2)),
        record$PARAMCD == "TRTCOMP" ~ as.numeric(round(runif(1, 0, 125), 0)),
        record$PARAMCD == "TNUMDOS" ~ as.numeric(round(runif(1, 0, 30), 0)),
        record$PARAMCD == "NUMCYC" ~ as.numeric(round(runif(1, 0, 24), 0)),
        record$PARAMCD == "MODEDD" ~ as.numeric(round(runif(1, 10, 50), 0)),
        record$PARAMCD == "FINDD" ~ as.numeric(round(runif(1, 10, 50), 0)),
        record$PARAMCD == "CUMDOSE" ~ as.numeric(round(runif(1, 100, 5000), 0)),
        record$PARAMCD == "DOSEDAYS" ~ as.numeric(round(runif(1, 1, 200), 0)),
        record$PARAMCD == "MEANDD" ~ as.numeric(round(runif(1, 5, 50), 1)),
        record$PARAMCD == "MEANDDI" ~ as.numeric(round(runif(1, 5, 45), 1))
      )

      # Add to records list
      records[[record_index]] <- record
      record_index <- record_index + 1
    }
  }

  # Combine all records
  gen <- dplyr::bind_rows(records)

  # Add additional columns
  gen <- dplyr::mutate(
    gen,
    AVALCAT1 = as.factor(dplyr::case_when(
      PARAMCD == "MODEDD" & AVAL <= 20 ~ "20 mg",
      PARAMCD == "MODEDD" & AVAL > 20 & AVAL <= 30 ~ "30 mg",
      PARAMCD == "MODEDD" & AVAL > 30 ~ "40 mg",
      PARAMCD == "FINDD" & AVAL <= 20 ~ "20 mg",
      PARAMCD == "FINDD" & AVAL > 20 & AVAL <= 30 ~ "30 mg",
      PARAMCD == "FINDD" & AVAL > 30 ~ "40 mg",
      PARAMCD == "TRTCOMP" & AVAL < 60 ~ "<60%",
      PARAMCD == "TRTCOMP" & AVAL >= 60 & AVAL < 80 ~ "60% to <80%",
      PARAMCD == "TRTCOMP" & AVAL >= 80 & AVAL < 100 ~ "80% to <100%",
      PARAMCD == "TRTCOMP" & AVAL >= 100 ~ ">=100%",
      PARAMCD == "NUMCYC" & AVAL >= 1 & AVAL < 10 ~ "1 to <10",
      PARAMCD == "NUMCYC" & AVAL >= 10 & AVAL < 20 ~ "10 to <20",
      PARAMCD == "NUMCYC" & AVAL >= 20 ~ ">=20",
      PARAMCD == "TNUMDOS" & AVAL >= 1 & AVAL < 10 ~ "1 to <10",
      PARAMCD == "TNUMDOS" & AVAL >= 10 & AVAL < 20 ~ "10 to <20",
      PARAMCD == "TNUMDOS" & AVAL >= 20 ~ ">=20",
      PARAMCD == "TRTDURM" & AVAL > 0 & AVAL < 3 ~ "0 to <3 months",
      PARAMCD == "TRTDURM" & AVAL >= 3 & AVAL < 6 ~ "3 to <6 months",
      PARAMCD == "TRTDURM" & AVAL >= 6 & AVAL < 9 ~ "6 to <9 months",
      PARAMCD == "TRTDURM" & AVAL >= 9 & AVAL < 12 ~ "9 to <12 months",
      PARAMCD == "TRTDURM" & AVAL >= 12 & AVAL < 15 ~ "12 to <15 months",
      PARAMCD == "TRTDURM" & AVAL >= 15 & AVAL < 18 ~ "15 to <18 months",
      PARAMCD == "TRTDURM" & AVAL >= 18 & AVAL < 21 ~ "18 to <21 months",
      PARAMCD == "TRTDURM" & AVAL >= 21 & AVAL < 24 ~ "21 to <24 months",
      PARAMCD == "TRTDURM" & AVAL >= 24 & AVAL < 27 ~ "24 to <27 months",
      PARAMCD == "TRTDURM" & AVAL >= 27 & AVAL < 30 ~ "27 to <30 months",
      PARAMCD == "TRTDURM" & AVAL >= 30 & AVAL < 33 ~ "30 to <33 months",
      PARAMCD == "TRTDURM" & AVAL >= 33 & AVAL < 36 ~ "33 to <36 months",
      PARAMCD == "TRTDURM" & AVAL >= 36 & AVAL < 39 ~ "36 to <39 months",
      .default = NA
    )),
    AVALCA1N = dplyr::case_when(
      PARAMCD == "MODEDD" & AVAL <= 20 ~ 1,
      PARAMCD == "MODEDD" & AVAL > 20 & AVAL <= 30 ~ 2,
      PARAMCD == "MODEDD" & AVAL > 30 ~ 3,
      PARAMCD == "FINDD" & AVAL <= 20 ~ 1,
      PARAMCD == "FINDD" & AVAL > 20 & AVAL <= 30 ~ 2,
      PARAMCD == "FINDD" & AVAL > 30 ~ 3,
      PARAMCD == "TRTCOMP" & AVAL < 60 ~ 1,
      PARAMCD == "TRTCOMP" & AVAL >= 60 & AVAL < 80 ~ 2,
      PARAMCD == "TRTCOMP" & AVAL >= 80 & AVAL < 100 ~ 3,
      PARAMCD == "TRTCOMP" & AVAL >= 100 ~ 4,
      PARAMCD == "NUMCYC" & AVAL >= 1 & AVAL < 10 ~ 1,
      PARAMCD == "NUMCYC" & AVAL >= 10 & AVAL < 20 ~ 2,
      PARAMCD == "NUMCYC" & AVAL >= 20 ~ 3,
      PARAMCD == "TNUMDOS" & AVAL >= 1 & AVAL < 10 ~ 1,
      PARAMCD == "TNUMDOS" & AVAL >= 10 & AVAL < 20 ~ 2,
      PARAMCD == "TNUMDOS" & AVAL >= 20 ~ 3,
      PARAMCD == "TRTDURM" & AVAL > 0 & AVAL < 3 ~ 1,
      PARAMCD == "TRTDURM" & AVAL >= 3 & AVAL < 6 ~ 2,
      PARAMCD == "TRTDURM" & AVAL >= 6 & AVAL < 9 ~ 3,
      PARAMCD == "TRTDURM" & AVAL >= 9 & AVAL < 12 ~ 4,
      PARAMCD == "TRTDURM" & AVAL >= 12 & AVAL < 15 ~ 5,
      PARAMCD == "TRTDURM" & AVAL >= 15 & AVAL < 18 ~ 6,
      PARAMCD == "TRTDURM" & AVAL >= 18 & AVAL < 21 ~ 7,
      PARAMCD == "TRTDURM" & AVAL >= 21 & AVAL < 24 ~ 8,
      PARAMCD == "TRTDURM" & AVAL >= 24 & AVAL < 27 ~ 9,
      PARAMCD == "TRTDURM" & AVAL >= 27 & AVAL < 30 ~ 10,
      PARAMCD == "TRTDURM" & AVAL >= 30 & AVAL < 33 ~ 11,
      PARAMCD == "TRTDURM" & AVAL >= 33 & AVAL < 36 ~ 12,
      PARAMCD == "TRTDURM" & AVAL >= 36 & AVAL < 39 ~ 13,
      .default = NA
    ),
    CRIT1 = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" ~ ">=1",
      PARAMCD == "TNUMDOS" ~ ">=1",
      PARAMCD == "TRTDURM" ~ "<1 month",
      .default = NA
    )),
    CRIT1FL = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" & AVAL >= 1 ~ "Y",
      PARAMCD == "NUMCYC" & AVAL < 1 ~ "N",
      PARAMCD == "TNUMDOS" & AVAL >= 1 ~ "Y",
      PARAMCD == "TNUMDOS" & AVAL < 1 ~ "N",
      PARAMCD == "TRTDURM" & AVAL < 1 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL >= 1 ~ "N",
      .default = NA
    )),
    CRIT2 = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" ~ ">=2",
      PARAMCD == "TNUMDOS" ~ ">=2",
      PARAMCD == "TRTDURM" ~ ">=1 month",
      .default = NA
    )),
    CRIT2FL = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" & AVAL >= 2 ~ "Y",
      PARAMCD == "NUMCYC" & AVAL < 2 ~ "N",
      PARAMCD == "TNUMDOS" & AVAL >= 2 ~ "Y",
      PARAMCD == "TNUMDOS" & AVAL < 2 ~ "N",
      PARAMCD == "TRTDURM" & AVAL >= 1 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL < 1 ~ "N",
      .default = NA
    )),
    CRIT3 = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" ~ ">=3",
      PARAMCD == "TNUMDOS" ~ ">=3",
      PARAMCD == "TRTDURM" ~ ">=2 months",
      .default = NA
    )),
    CRIT3FL = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" & AVAL >= 3 ~ "Y",
      PARAMCD == "NUMCYC" & AVAL < 3 ~ "N",
      PARAMCD == "TNUMDOS" & AVAL >= 3 ~ "Y",
      PARAMCD == "TNUMDOS" & AVAL < 3 ~ "N",
      PARAMCD == "TRTDURM" & AVAL >= 2 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL < 2 ~ "N",
      .default = NA
    )),
    CRIT4 = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" ~ ">=5",
      PARAMCD == "TNUMDOS" ~ ">=5",
      PARAMCD == "TRTDURM" ~ ">=3 months",
      .default = NA
    )),
    CRIT4FL = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" & AVAL >= 5 ~ "Y",
      PARAMCD == "NUMCYC" & AVAL < 5 ~ "N",
      PARAMCD == "TNUMDOS" & AVAL >= 5 ~ "Y",
      PARAMCD == "TNUMDOS" & AVAL < 5 ~ "N",
      PARAMCD == "TRTDURM" & AVAL >= 3 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL < 3 ~ "N",
      .default = NA
    )),
    CRIT5 = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" ~ ">=10",
      PARAMCD == "TNUMDOS" ~ ">=10",
      PARAMCD == "TRTDURM" ~ ">=6 months",
      .default = NA
    )),
    CRIT5FL = as.factor(dplyr::case_when(
      PARAMCD == "NUMCYC" & AVAL >= 10 ~ "Y",
      PARAMCD == "NUMCYC" & AVAL < 10 ~ "N",
      PARAMCD == "TNUMDOS" & AVAL >= 10 ~ "Y",
      PARAMCD == "TNUMDOS" & AVAL < 10 ~ "N",
      PARAMCD == "TRTDURM" & AVAL >= 6 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL < 6 ~ "N",
      .default = NA
    )),
    CRIT6 = as.factor(dplyr::case_when(
      PARAMCD == "TRTDURM" ~ ">=12 months",
      .default = NA
    )),
    CRIT6FL = as.factor(dplyr::case_when(
      PARAMCD == "TRTDURM" & AVAL >= 12 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL < 12 ~ "N",
      .default = NA
    )),
    CRIT7 = as.factor(dplyr::case_when(
      PARAMCD == "TRTDURM" ~ ">=24 months",
      .default = NA
    )),
    CRIT7FL = as.factor(dplyr::case_when(
      PARAMCD == "TRTDURM" & AVAL >= 24 ~ "Y",
      PARAMCD == "TRTDURM" & AVAL < 24 ~ "N",
      .default = NA
    )),
    AVISIT = "Overall",
    AVISITN = 1
  )

  # make sure the levels matchs
  gen$PARAMCD <- factor(
    gen$PARAMCD,
    levels = c(
      "CUMDOSE",
      "DINTEN",
      "DOSEDAYS",
      "FINDD",
      "MEANDD",
      "MEANDDI",
      "MODEDD",
      "NUMCYC",
      "RDINTE",
      "TNUMDOS",
      "TRTCOMP",
      "TRTDURM",
      "TRTDURY"
    )
  )

  # Select only needed columns
  gen <- dplyr::select(
    gen,
    USUBJID,
    PARAMCD,
    PARAM,
    AVAL,
    AVALCAT1,
    AVALCA1N,
    CRIT1,
    CRIT1FL,
    CRIT2,
    CRIT2FL,
    CRIT3,
    CRIT3FL,
    CRIT4,
    CRIT4FL,
    CRIT5,
    CRIT5FL,
    CRIT6,
    CRIT6FL,
    CRIT7,
    CRIT7FL,
    AVISIT,
    AVISITN
  )

  # Ensure one record per subject per parameter
  gen <- dplyr::group_by(gen, USUBJID, PARAMCD)
  gen <- dplyr::slice(gen, 1)
  gen <- dplyr::ungroup(gen)

  source(file.path("data-raw", "adsl.R"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), "USUBJID")

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c("TRT01A", "SAFFL", "STUDYID")

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl %>%
    select(USUBJID, all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    message("Dropping shared vars from raw: ", paste(shared, collapse = ", "))
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID")

  # Define additional labels
  additional_labels <- list(
    AVISIT = "Analysis Visit",
    AVISITN = "Analysis Visit (N)",
    USUBJID = "Unique Subject Identifier",
    PARAM = "Parameter",
    AVAL = "Analysis Value",
    AVALCAT1 = "Analysis Value Category 1",
    AVALCA1N = "Analysis Value Category 1 (N)",
    CRIT1 = "Analysis Criterion 1",
    CRIT1FL = "Criterion 1 Evaluation Result Flag",
    CRIT2 = "Analysis Criterion 2",
    CRIT2FL = "Criterion 2 Evaluation Result Flag",
    CRIT3 = "Analysis Criterion 3",
    CRIT3FL = "Criterion 3 Evaluation Result Flag",
    CRIT4 = "Analysis Criterion 4",
    CRIT4FL = "Criterion 4 Evaluation Result Flag",
    CRIT5 = "Analysis Criterion 5",
    CRIT5FL = "Criterion 5 Evaluation Result Flag",
    CRIT6 = "Analysis Criterion 6",
    CRIT6FL = "Criterion 6 Evaluation Result Flag",
    CRIT7 = "Analysis Criterion 7",
    CRIT7FL = "Criterion 7 Evaluation Result Flag",
    PARAMCD = "Parameter Code"
  )

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Apply labels directly
  for (var_name in names(additional_labels)) {
    if (var_name %in% names(gen)) {
      attr(gen[[var_name]], "label") <- additional_labels[[var_name]]
    }
  }

  return(gen)
}

# Generate and save the dataset
adexsum <- gen_adexsum()
