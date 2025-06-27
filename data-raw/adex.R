# Generate ADEX dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(formatters)
library(forcats)
library(stringr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADEX dataset
gen_adex <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  # Get source data
  raw <- pharmaverseadam::adex

  raw <- dplyr::filter(raw, PARAMCD == "DOSE")

  gen <- dplyr::mutate(
    raw,
    EXTRT = forcats::fct_recode(
      EXTRT,
      "APALUTAMIDE" = "XANOMELINE"
    ),
    EXDOSE = dplyr::case_when(
      EXDOSE == 54 ~ 7.5,
      EXDOSE == 81 ~ 15,
      .default = EXDOSE
    ),
    TRT01P = forcats::fct_recode(
      TRT01P,
      "Apalutamide" = "Xanomeline High Dose",
      "Apalutamide Subgroup" = "Xanomeline Low Dose"
    ),
    ATRT = as.factor(dplyr::case_when(
      TRT01P == "Apalutamide" ~ "APALUTAMIDE",
      TRT01P == "Apalutamide Subgroup" ~ "APALUTAMIDE",
      TRT01P == "Placebo" ~ "PLACEBO"
    )),
    ARMCD = as.factor(dplyr::case_when(
      ARMCD == "Xan_Hi" ~ "Apa",
      ARMCD == "Xan_Lo" ~ "Apa_Sub",
      .default = ARMCD
    )),
    ARM = as.factor(dplyr::case_when(
      ARM == "Xanomeline High Dose" ~ "Apalutamide",
      ARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
      .default = ARM
    )),
    ACTARMCD = as.factor(dplyr::case_when(
      ACTARMCD == "Xan_Hi" ~ "Apa",
      ACTARMCD == "Xan_Lo" ~ "Apa_Sub",
      .default = ACTARMCD
    )),
    ACTARM = as.factor(dplyr::case_when(
      ACTARM == "Xanomeline High Dose" ~ "Apalutamide",
      ACTARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
      .default = ACTARM
    )),
    DAEXPDTC = as.Date(sample(
      c("2013-09-10", "2013-12-15", "2014-02-05", "2014-03-20"),
      dplyr::n(),
      replace = TRUE
    )),
    EXLOT = case_when(
      DAEXPDTC == as.Date("2013-09-10") ~ "XXX-YYY-ZZZ-004",
      DAEXPDTC == as.Date("2013-12-15") ~ "XXX-YYY-ZZZ-002",
      DAEXPDTC == as.Date("2014-02-05") ~ "XXX-YYY-ZZZ-003",
      DAEXPDTC == as.Date("2014-03-20") ~ "XXX-YYY-ZZZ-005",
      TRUE ~ "XXX-YYY-ZZZ-006" # Default value for all other dates
    ),
    ADOSE = case_when(
      is.na(DOSEO) ~ NA_real_,
      DOSEO == 0 ~ 60,
      DOSEO < 60 ~ 60,
      DOSEO >= 60 & DOSEO < 180 ~ 120,
      DOSEO >= 180 & DOSEO < 240 ~ 180,
      DOSEO >= 240 & DOSEO < 300 ~ 240,
      DOSEO >= 300 & DOSEO < 480 ~ 300,
      DOSEO >= 480 & DOSEO < 720 ~ 480,
      DOSEO >= 720 & DOSEO < 1000 ~ 720,
      DOSEO >= 1000 & DOSEO < 5000 ~ 960,
      DOSEO >= 5000 & DOSEO < 10000 ~ 1800,
      DOSEO >= 10000 ~ 3600
    ),
    TRT01P = droplevels(dplyr::case_when(
      TRT01P == "Screen Failure" ~ NA,
      .default = TRT01P
    )),
    TRT01PN = dplyr::case_when(
      TRT01P == "Apalutamide" ~ 1,
      TRT01P == "Apalutamide Subgroup" ~ 2,
      TRT01P == "Placebo" ~ 3
    ),
    TRT01P = forcats::fct_reorder(TRT01P, TRT01PN, .na_rm = TRUE),
    TRT01A = forcats::fct_recode(
      TRT01A,
      "Apalutamide" = "Xanomeline High Dose",
      "Apalutamide Subgroup" = "Xanomeline Low Dose"
    ),
    TRT01A = droplevels(dplyr::case_when(
      TRT01A == "Screen Failure" ~ NA,
      .default = TRT01A
    )),
    AVISITN = case_when(
      VISIT == "Baseline" ~ 1,
      VISIT == "Week 2" ~ 2,
      VISIT == "Week 4" ~ 3,
      VISIT == "Week 6" ~ 4,
      VISIT == "Week 8" ~ 5,
      VISIT == "Week 12" ~ 6,
      VISIT == "Week 16" ~ 7,
      VISIT == "Week 20" ~ 8,
      VISIT == "Week 24" ~ 9,
      VISIT == "Week 26" ~ 10,
    ),
    AVISIT = fct_reorder(
      as.factor(case_when(
        VISIT == "Baseline" ~ "Screening",
        VISIT == "Week 2" ~ "Cycle 02",
        VISIT == "Week 4" ~ "Cycle 03",
        VISIT == "Week 6" ~ "Cycle 04",
        VISIT == "Week 8" ~ "Cycle 05",
        VISIT == "Week 12" ~ "Cycle 06",
        VISIT == "Week 16" ~ "Cycle 07",
        VISIT == "Week 20" ~ "Cycle 08",
        VISIT == "Week 24" ~ "Cycle 09",
        VISIT == "Week 26" ~ "End Of Treatment",
        TRUE ~ as.character(VISIT) # Other values remain unchanged
      )),
      AVISITN,
      .na_rm = FALSE
    ),
    TRT01AN = dplyr::case_when(
      TRT01A == "Apalutamide" ~ 1,
      TRT01A == "Apalutamide Subgroup" ~ 2,
      TRT01A == "Placebo" ~ 3
    ),
    TRT01A = forcats::fct_reorder(TRT01A, TRT01AN, .na_rm = TRUE),
    AGEGR1 = as.factor(dplyr::case_when(
      AGE >= 18 & AGE < 65 ~ ">=18 to <65",
      AGE >= 65 & AGE < 75 ~ ">=65 to <75",
      AGE >= 75 ~ ">=75"
    )),
    AOCCUR = as.factor(sample(c("N", "Y"), dplyr::n(), replace = TRUE)),
    SEX = as.factor(dplyr::case_when(
      SEX == "F" ~ "Female",
      SEX == "M" ~ "Male"
    )),
    COUNTRY = as.factor("United States of America"),
    RACE = as.factor(dplyr::case_when(
      RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~
        "American Indian or Alaska Native",
      RACE == "ASIAN" ~ "Asian",
      RACE == "BLACK OR AFRICAN AMERICAN" ~ "Black or African American",
      RACE == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~
        "Native Hawaiian or other Pacific Islander",
      RACE == "WHITE" ~ "White",
      RACE == "MULTIPLE" ~ "Multiple",
      RACE == "NOT REPORTED" ~ "Not reported",
      RACE == "UNKNOWN" ~ "Unknown",
      RACE == "OTHER" ~ "Other",
    )),
    RACE_DECODE = as.factor(dplyr::case_when(
      RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~
        "American Indian or Alaska Native",
      RACE == "ASIAN" ~ "Asian",
      RACE == "BLACK OR AFRICAN AMERICAN" ~ "Black or African American",
      RACE == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~
        "Native Hawaiian or other Pacific Islander",
      RACE == "WHITE" ~ "White",
      RACE == "MULTIPLE" ~ "Multiple",
      RACE == "NOT REPORTED" ~ "Not reported",
      RACE == "UNKNOWN" ~ "Unknown",
      RACE == "OTHER" ~ "Other"
    )),
    RACEGR1 = as.factor(RACEGR1),
    ETHNIC = as.factor(dplyr::case_when(
      ETHNIC == "HISPANIC OR LATINO" ~ "Hispanic or Latino",
      ETHNIC == "NOT HISPANIC OR LATINO" ~ "Not Hispanic or Latino",
      ETHNIC == "NOT REPORTED" ~ "Not reported",
      ETHNIC == "UNKNOWN" ~ "Unknown",
    )),
    ACAT1 = as.factor(case_when(
      is.na(EXADJ) ~ "Dose not adjusted",
      EXADJ == "ADVERSE EVENT" ~ "Dose adjusted",
      EXADJ == "MEDICAL ERROR" ~ "Dose adjusted",
      TRUE ~ "Other" # default case if needed
    )),
    AREASOC = as.factor(case_when(
      is.na(EXADJ) ~ "",
      EXADJ == "ADVERSE EVENT" ~ "Adverse Event",
      EXADJ == "MEDICATION ERROR" ~ "Other",
      TRUE ~ "" # default case
    )),
    AREASOO = as.factor(case_when(
      AREASOC == "Other" ~ "Other reason",
      TRUE ~ "" # default case
    )),
    AADJ = AREASOC,
    AADJPOTH = dplyr::case_when(
      AADJ == "Other" ~ "************",
      .default = "Reason prior to infusion"
    ),
    AADJP = AADJ,
    AACTDU = as.factor(sample(
      c(
        "INFUSION INTERRUPTED",
        "",
        "INFUSION RATE INCREASED",
        "INFUSION CONTINUED AT SAME RATE",
        "INFUSION ABORTED",
        "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AACTDU1 = as.factor(case_when(
      AACTDU == "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE" ~
        "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE",
      TRUE ~ "" # default case
    )),
    AACTDU2 = case_when(
      AACTDU == "INFUSION ABORTED" ~ "INFUSION ABORTED",
      TRUE ~ "" # default case
    ),
    AACTDU3 = case_when(
      AACTDU == "INFUSION INTERRUPTED" ~ "INFUSION INTERRUPTED",
      TRUE ~ "" # default case
    ),
    AACTDU4 = case_when(
      AACTDU == "INFUSION RATE DECREASED" ~ "INFUSION RATE DECREASED",
      TRUE ~ "" # default case
    ),
    AACTDU5 = case_when(
      AACTDU == "INFUSION RATE INCREASED" ~ "INFUSION RATE INCREASED",
      TRUE ~ "" # default case
    ),
    AADJOTH = ifelse(
      AREASOC == "Other",
      sample(
        c(
          "CYSTOSCOPI WITH LESION. RTU PENDING",
          "MULTIPLE TUMORS IDENTIFIED BY CYSTOSCOPY ON WEEK 12.",
          "DURING WEEK 12, THE CYSTOSCOPY WAS PERFORMED AND SUSPICION OF PROGRESSION WAS OBSERVED.",
          "BLADDER MAPPING PERFORMED ON 20 NOV 2023",
          "PROGRESSION",
          "DRUG WAS NEVER ADMINISTERED BECAUSE THE SUBJECT DROPPED OUT PRIOR TO FIRST DOSE",
          "PERSONAL REASONS",
          "TAR200 WAS HALF BLOCKED BY THE EXIT PORT OF THE UPC",
          "PQC",
          "THE DEVICE COULD NOT PUSH TAR-200 OUT.",
          "SHEET COULD NOT INTO THE BLADDER",
          "PT HAD HIS 24W TURB ON 22TH OF DECEMBER SO TAR WAS REMOVED AT THAT TIME",
          "THE DOSE DELAYED DUE TO TURB PERFORMED 10-NOV-2023",
          "TAR-200 DID NOT TAKE THE PRETZEL SHAPE INSIDE THE BLADDER.",
          "PATIENT'S VACATION",
          "TAR200 WAS BLOCKED IN THE MIDDLE OF THE UPC HOLE",
          "PQC WAS REPORTED.",
          "PQC WAS IDENTIFIED",
          "URINARY PLACEMENT CATHETER COULD NOT INTO THE BLADDER.",
          "INSERTION FAILED",
          "SUSPICION OF DISEASE PROGRESSION",
          "DUE TO PERSONAL SCHEDULE",
          "URINE CULTURE NOT DISPONIBLE"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      "Reason during infusion"
    ),
    ACAT2 = as.factor(sample(
      c("Dose not administered", "Dose administered"),
      dplyr::n(),
      replace = TRUE
    )),
    AACTPR = as.factor(sample(
      c(
        "INFUSION RATE DECREASED COMPARED TO PRIOR INFUSION",
        "",
        "INFUSION SKIPPED (AND NOT MADE UP)",
        "DOSE RE-ESCALATED",
        "INFUSION DELAYED WITHIN THE CYCLE",
        "DOSE REDUCED COMPARED TO PRIOR INFUSION",
        "SAME DOSE AS PRIOR INFUSION",
        "STUDY DRUG PERMANENTLY DISCONTINUED"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AACTPR_DECODE = stringr::str_to_sentence(AACTPR),
    ASCHDOSE = EXDOSE,
    ASCHDOSU = EXDOSU,
    ADOSFRM = EXDOSFRM,
    ADOSU = EXDOSU,
    ADOSFRQ = EXDOSFRQ,
    AROUTE = EXROUTE,
    ATVINF = ADOSE,
    ATVINFU = ADOSU,
    AINFRAT = ADOSE,
    AINFRAU = ADOSU
  )

  # Additional labels for all relevant variables
  additional_labels <- list(
    EXTRT = "Planned Treatment",
    EXDOSE = "Adjusted Dose",
    ACAT1 = "Analysis Category 1",
    AADJOTH = "Other Anal Reason for Dose Adjustment",
    ACAT2 = "Analysis Category 2",
    AACTPR = "Action Taken Prior to Infusion Start",
    AREASOC = "Analysis Reason for Occur Value",
    AADJ = "Analysis Reason for Dose Adjustment",
    AADJP = "Analysis Reason for Dose Adjustment Prior",
    AREASOO = "Other Analysis Reason for Occur Value",
    AACTDU = "Analysis Action Taken During Study Trt",
    AACTDU1 = "Act Takn Dur Infus-Full Dose Admined",
    AACTDU2 = "Act Takn Dur Infus-Infusion Aborted",
    AACTDU3 = "Act Takn Dur Infus-Infusion Interrupted",
    AACTDU4 = "Act Takn Dur Infus-Infusion Rate Decrsed",
    AACTDU5 = "Act Takn Dur Infus-Infusion Rate Incrsed",
    TRT01P = "Planned Treatment for Period 01",
    TRT01A = "Actual Treatment for Period 01",
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    ARMCD = "Treatment Category Code",
    ARM = "Treatment Group",
    ACTARMCD = "Actual Arm Code",
    ACTARM = "Actual Treatment Group",
    DAEXPDTC = "Date of Exposure",
    EXLOT = "Lot Number",
    ADOSE = "Analysis Dose",
    AVISITN = "Visit Number",
    AVISIT = "Visit Label",
    ADOSE = "Analysis Dose",
    AGEGR1 = "Age Group",
    AOCCUR = "Analysis Occurrence",
    SEX = "Sex",
    COUNTRY = "Country",
    RACE = "Race",
    ETHNIC = "Ethnicity",
    RACE_DECODE = "Race",
    ATRT = "Analysis name of Treatment",
    ASCHDOSE = "Analysis Scheduled Dose",
    ASCHDOSU = "Analysis Scheduled Dose Units",
    ADOSFRM = "Analysis Dose Form",
    ADOSU = "Analysis Dose Units",
    ADOSFRQ = "Analysis Dosing Frequency per Interval",
    AROUTE = "Analysis Route of Administration",
    AADJPOTH = "Other Anal Reason for Dose Adjust Prior",
    AACTPR_DECODE = "Action Taken Prior to Infusion Start",
    ATVINF = "Analysis Total Volume Infused",
    ATVINFU = "Analysis Total Volume Infused Units",
    AINFRAT = "Analysis Infusion Rate",
    AINFRAU = "Analysis Infusion Rate Unit"
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
adex <- gen_adex()
