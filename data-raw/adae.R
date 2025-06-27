# Generate ADAE dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(formatters)
library(forcats)
library(admiral)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADAE dataset
gen_adae <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  # Get source data
  raw <- pharmaverseadam::adae

  gen <- dplyr::mutate(
    raw,
    AETOXGR = as.factor(sample(seq(0, 5), dplyr::n(), replace = TRUE)),
    AETOXGRN = as.numeric(as.character(AETOXGR)),
    AEACN = as.factor(sample(
      c(
        "DOSE NOT CHANGED",
        "NOT APPLICABLE",
        "DRUG WITHDRAWN",
        "DOSE REDUCED",
        "DRUG INTERRUPTED"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AEACN_DECODE = dplyr::case_when(
      AEACN == "DOSE NOT CHANGED" ~ "Dose Not Changed",
      AEACN == "NOT APPLICABLE" ~ "Not Applicable",
      AEACN == "DRUG WITHDRAWN" ~ "Drug Withdrawn",
      AEACN == "DOSE REDUCED" ~ "Dose Reduced",
      AEACN == "DRUG INTERRUPTED" ~ "Drug Interrupted"
    ),
    AESEV = dplyr::case_when(
      AESEV == "MILD" ~ "Mild",
      AESEV == "MODERATE" ~ "Moderate",
      AESEV == "SEVERE" ~ "Severe"
    ),
    DOSEDY = as.numeric(37),
    DOSEU = as.factor("mg"),
    DOSEON = as.numeric(10),
    AECONTRT = as.factor(sample(c("N", "Y", "U"), dplyr::n(), replace = TRUE)),
    CQ01NAM = as.factor(sample(c("Seizure", NA), dplyr::n(), replace = TRUE)),
    CQ02NAM = as.factor(sample(c("Skin rash", NA), dplyr::n(), replace = TRUE)),
    CQ03NAM = as.factor(sample(
      c("Hypothyroidism", NA),
      dplyr::n(),
      replace = TRUE
    )),
    AESMIE = "Y",
    AESER = as.factor(sample(c("N", "Y"), dplyr::n(), replace = TRUE)),
    AESER_DECODE = dplyr::case_when(
      AESER == "Y" ~ "Yes",
      AESER == "N" ~ "No",
      .default = NA
    ),
    AEREL = as.factor(dplyr::case_when(
      AEREL == "PROBABLE" ~ "PROBABLE",
      AEREL == "REMOTE" ~ "RELATED",
      AEREL == "POSSIBLE" ~ "POSSIBLE",
      AEREL == "NONE" ~ "NOT RELATED",
      is.na(AEREL) ~ NA_character_
    )),
    AEREL_DECODE = as.factor(dplyr::case_when(
      AEREL == "PROBABLE" ~ "Probable",
      AEREL == "REMOTE" ~ "Remote",
      AEREL == "POSSIBLE" ~ "Possible",
      AEREL == "NONE" ~ "None",
      is.na(AEREL) ~ "Not applicable"
    )),
    AEOUT_DECODE = as.factor(dplyr::case_when(
      AEOUT == "NOT RECOVERED/NOT RESOLVED" ~ "Not Recovered/Not Resolved",
      AEOUT == "RECOVERED/RESOLVED" ~ "Recovered/Resolved",
      AEOUT == "FATAL" ~ "Fatal",
      TRUE ~ "Other"
    )),
    AEBODSYS = forcats::fct_relabel(AEBODSYS, stringr::str_to_sentence) # Convert AEBODSYS levels to sentence
  ) %>%
    derive_var_extreme_flag(
      new_var = AOCCFL,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(STUDYID, USUBJID, ASTDY, AESEQ),
      mode = "first"
    ) %>%
    derive_var_extreme_flag(
      new_var = AOCCPFL,
      by_vars = exprs(STUDYID, USUBJID, AEDECOD),
      order = exprs(STUDYID, USUBJID, AEDECOD, ASTDY, AESEQ),
      mode = "first"
    ) %>%
    derive_var_extreme_flag(
      new_var = AOCCSFL,
      by_vars = exprs(STUDYID, USUBJID, AEBODSYS),
      order = exprs(STUDYID, USUBJID, AEDECOD, ASTDY, AESEQ),
      mode = "first"
    )

  source(file.path("data-raw", "adsl.R"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), "USUBJID")

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c(
    "TRT01A",
    "SAFFL",
    "AGE",
    "SEX",
    "RACE",
    "RACE_DECODE",
    "STUDYID",
    "AGEGR1",
    "TRTEDY"
  )

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl %>%
    select(USUBJID, all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    message("Dropping shared vars from raw: ", paste(shared, collapse = ", "))
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID") %>%
    mutate(
      months = (TRTEDY + 30) / 30.4375,
      ACAT1 = case_when(
        months <= 3 ~ "Within 3 months",
        months > 3 & months <= 6 ~ "4 to 6 months",
        months > 6 & months <= 9 ~ "7 to 9 months",
        months > 9 & months <= 12 ~ "10 to 12 months",
        months > 12 ~ "Beyond 13 months",
        .default = NA_character_
      ),
      ACAT1 = factor(ACAT1, levels = c("Within 3 months", "4 to 6 months", "7 to 9 months", "10 to 12 months", "Beyond 13 months"))
    )

  # Add labels
  additional_labels <- list(
    SAFFL = "Safety Population Flag",
    AESER = "Serious Event",
    AESER_DECODE = "Serious Event",
    ACAT1 = "Analysis Category 1",
    AETOXGR = "Standard Toxicity Grade",
    AETOXGRN = "Standard Toxicity Grade (N)",
    AEACN_DECODE = "Action Taken with Study Treatment",
    DOSEDY = "Day of Study Drug",
    DOSEU = "Treatment Dose Units",
    DOSEON = "Treatment Dose at Record Start",
    AEREL_DECODE = "Causality",
    AEOUT_DECODE = "Outcome of Adverse Event",
    AOCCFL = "1st Occurance within Subject Flag",
    AOCCPFL = "1st Occurance within Preferred Term Flag",
    AOCCSFL = "1st Occurrence of SOC Flag",
    CQ01NAM = "Customized Query 01 Name",
    CQ02NAM = "Customized Query 02 Name",
    CQ03NAM = "Customized Query 03 Name",
    AECONTRT = "Concomitant or Additional Trtmnt Given",
    AENDTF_DECODE = "Analysis End Date Imputation Flag",
    AESMIE = "Other Medically Important Serious Event",
    AESMIE_DECODE = "Other Medically Important Serious Event"
  )

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels
  gen <- restore_labels(
    df = gen,
    orig_df = raw,
    additional_labels = additional_labels
  )

  return(gen)
}

# Generate the dataset
adae <- gen_adae()
