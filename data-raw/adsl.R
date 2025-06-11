# Generate ADSL dataset

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(formatters)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADSL dataset
gen_adsl <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get source data
  raw <- pharmaverseadam::adsl

  gen <- raw
  gen$TRT01P <- forcats::fct_recode(
    gen$TRT01P,
    "Apalutamide" = "Xanomeline High Dose",
    "Apalutamide Subgroup" = "Xanomeline Low Dose"
  )
  gen$ARMCD <- as.factor(dplyr::case_when(
    gen$ARMCD == "Xan_Hi" ~ "Apa",
    gen$ARMCD == "Xan_Lo" ~ "Apa_Sub",
    .default = gen$ARMCD
  ))
  gen$ARM <- as.factor(dplyr::case_when(
    gen$ARM == "Xanomeline High Dose" ~ "Apalutamide",
    gen$ARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
    .default = gen$ARM
  ))
  gen$ACTARMCD <- as.factor(dplyr::case_when(
    gen$ACTARMCD == "Xan_Hi" ~ "Apa",
    gen$ACTARMCD == "Xan_Lo" ~ "Apa_Sub",
    .default = gen$ACTARMCD
  ))
  gen$ACTARM <- as.factor(dplyr::case_when(
    gen$ACTARM == "Xanomeline High Dose" ~ "Apalutamide",
    gen$ACTARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
    .default = gen$ACTARM
  ))
  gen$TRT01P <- droplevels(dplyr::case_when(
    gen$TRT01P == "Screen Failure" ~ NA,
    .default = gen$TRT01P
  ))
  gen$TRT01PN <- dplyr::case_when(
    gen$TRT01P == "Apalutamide" ~ 1,
    gen$TRT01P == "Apalutamide Subgroup" ~ 2,
    gen$TRT01P == "Placebo" ~ 3
  )
  gen$TRT01P <- forcats::fct_reorder(gen$TRT01P, gen$TRT01PN, .na_rm = TRUE)
  gen$TRT01A <- forcats::fct_recode(
    gen$TRT01A,
    "Apalutamide" = "Xanomeline High Dose",
    "Apalutamide Subgroup" = "Xanomeline Low Dose"
  )
  gen$TRT01A <- droplevels(dplyr::case_when(
    gen$TRT01A == "Screen Failure" ~ NA,
    .default = gen$TRT01A
  ))
  gen$TRT01AN <- dplyr::case_when(
    gen$TRT01A == "Apalutamide" ~ 1,
    gen$TRT01A == "Apalutamide Subgroup" ~ 2,
    gen$TRT01A == "Placebo" ~ 3
  )
  gen$TRT01A <- forcats::fct_reorder(gen$TRT01A, gen$TRT01AN, .na_rm = TRUE)
  gen$AGEGR1 <- as.factor(dplyr::case_when(
    gen$AGE >= 18 & gen$AGE < 65 ~ ">=18 to <65",
    gen$AGE >= 65 & gen$AGE < 75 ~ ">=65 to <75",
    gen$AGE >= 75 ~ ">=75"
  ))
  gen$AGEGR1N <- dplyr::case_when(
    gen$AGEGR1 == ">=18 to <65" ~ 1,
    gen$AGEGR1 == ">=65 to <75" ~ 2,
    gen$AGEGR1 == ">=75" ~ 3
  )
  gen$SEX_DECODE <- as.factor(dplyr::case_when(
    gen$SEX == "F" ~ "Female",
    gen$SEX == "M" ~ "Male"
  ))
  gen$WEIGHTBL <- as.numeric(sample(seq(0, 150), nrow(gen), replace = TRUE))
  gen$WGTGR1N <- dplyr::case_when(
    gen$WEIGHTBL < 30 ~ 1,
    gen$WEIGHTBL >= 30 & gen$WEIGHTBL < 60 ~ 2,
    gen$WEIGHTBL >= 60 & gen$WEIGHTBL < 90 ~ 3,
    gen$WEIGHTBL >= 90 ~ 4
  )
  gen$WGTGR1 <- forcats::fct_reorder(
    as.factor(dplyr::case_when(
      gen$WEIGHTBL < 30 ~ "<30",
      gen$WEIGHTBL >= 30 & gen$WEIGHTBL < 60 ~ "30 to <60",
      gen$WEIGHTBL >= 60 & gen$WEIGHTBL < 90 ~ "60 to <90",
      gen$WEIGHTBL >= 90 ~ ">=90"
    )),
    gen$WGTGR1N
  )
  gen$HEIGHTBL <- as.numeric(sample(seq(0, 150), nrow(gen), replace = TRUE))
  gen$BSABL <- as.numeric(sample(seq(0, 150), nrow(gen), replace = TRUE))
  gen$BMIBL <- as.numeric(sample(seq(0, 150), nrow(gen), replace = TRUE))
  gen$BMIBLG1N <- dplyr::case_when(
    gen$BMIBL < 18.5 ~ 1,
    gen$BMIBL >= 18.5 & gen$BMIBL < 25 ~ 2,
    gen$BMIBL >= 25 & gen$BMIBL < 30 ~ 3,
    gen$BMIBL >= 30 ~ 4
  )
  gen$BMIBLG1 <- forcats::fct_reorder(
    as.factor(dplyr::case_when(
      gen$BMIBL < 18.5 ~ "Underweight <18.5",
      gen$BMIBL >= 18.5 & gen$BMIBL < 25 ~ "Normal 18.5 to <25",
      gen$BMIBL >= 25 & gen$BMIBL < 30 ~ "Overweight 25 to <30",
      gen$BMIBL >= 30 ~ "Obese >=30"
    )),
    gen$BMIBLG1N
  )
  gen$SEX <- as.factor(gen$SEX)
  gen$RACE <- as.factor(gen$RACE)
  gen$ETHNIC <- as.factor(gen$ETHNIC)
  gen$COUNTRY_DECODE <- as.factor("United States of America")
  gen$RACE_DECODE <- as.factor(dplyr::case_when(
    gen$RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~
      "American Indian or Alaska Native",
    gen$RACE == "ASIAN" ~ "Asian",
    gen$RACE == "BLACK OR AFRICAN AMERICAN" ~ "Black or African American",
    gen$RACE == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~
      "Native Hawaiian or other Pacific Islander",
    gen$RACE == "WHITE" ~ "White",
    gen$RACE == "MULTIPLE" ~ "Multiple",
    gen$RACE == "NOT REPORTED" ~ "Not reported",
    gen$RACE == "UNKNOWN" ~ "Unknown",
    gen$RACE == "OTHER" ~ "Other"
  ))
  gen$RACEGR1 <- as.factor(gen$RACEGR1)
  gen$RFICDTC <- gen$DMDTC
  gen$RFICDT <- as.Date(gen$DMDTC)
  gen$ETHNIC_DECODE <- as.factor(dplyr::case_when(
    gen$ETHNIC == "HISPANIC OR LATINO" ~ "Hispanic or Latino",
    gen$ETHNIC == "NOT HISPANIC OR LATINO" ~ "Not Hispanic or Latino",
    gen$ETHNIC == "NOT REPORTED" ~ "Not reported",
    gen$ETHNIC == "UNKNOWN" ~ "Unknown"
  ))
  gen$STRAT1R <- as.factor("Stratification Factor 1")
  gen$STRAT2R <- as.factor("Stratification Factor 2")
  gen$RANUM <- as.factor("1000001")
  gen$RANDDTM <- as.POSIXct(
    paste0(gen$RANDDT, " 11:59"),
    format = "%Y-%m-%d %H:%M"
  )
  gen$EOTSTT <- as.factor(gen$EOSSTT)
  gen$DCTREAS <- as.factor(dplyr::case_when(
    gen$EOTSTT == "DISCONTINUED" ~ "Other",
    .default = NA
  ))
  gen$LTVISIT <- as.factor("Last Treatment Visit")
  gen$DCTREASP <- dplyr::case_when(
    gen$DCTREAS == "Other" ~ "specify text",
    .default = NA
  )
  gen$DCTDT <- dplyr::case_when(
    gen$EOTSTT == "DISCONTINUED" ~ gen$EOSDT,
    .default = NA
  )
  gen$DCSREAS <- as.factor(dplyr::case_when(
    gen$EOSSTT == "DISCONTINUED" ~ "Other",
    .default = NA
  ))
  gen$DCSREASP <- dplyr::case_when(
    gen$DCSREAS == "Other" ~ "specify text",
    .default = NA
  )
  gen$LSVISIT <- as.factor("Last Study Visit")

  # Generate TRTEDY to ensure all ACAT1 groups are represented correctly
  n <- nrow(gen)
  # Calculate max TRTEDY for each group so that (TRTEDY + 30)/30.4375 falls in the right range
  acat1_breaks <- list(
    rep(50, ceiling(n * 0.2)), # Within 3 months (TRTEDY=50)
    rep(100, ceiling(n * 0.2)), # 4 to 6 months (TRTEDY=100)
    rep(200, ceiling(n * 0.2)), # 7 to 9 months (TRTEDY=200)
    rep(300, ceiling(n * 0.2)), # 10 to 12 months (TRTEDY=300)
    rep(400, n - 4 * ceiling(n * 0.2)) # Beyond 13 months (TRTEDY=400)
  )
  trtedy_vals <- unlist(acat1_breaks)
  trtedy_vals <- sample(trtedy_vals, n)
  gen$TRTEDY <- trtedy_vals

  gen$SCRNFL <- as.factor("Y")
  gen$SCRFFL <- as.factor(dplyr::case_when(
    !is.na(gen$SCRFDT) ~ "Y",
    .default = "N"
  ))
  gen$DCSCREEN <- as.factor(dplyr::case_when(
    !is.na(gen$SCRFDT) ~ "Failure to meet eligibility criteria",
    .default = NA
  ))
  gen$SAFFL <- dplyr::case_when(
    is.na(gen$SAFFL) ~ "N",
    .default = gen$SAFFL
  )
  gen$ENRLFL <- as.factor(dplyr::case_when(
    gen$SAFFL == "Y" ~ "Y",
    .default = "N"
  ))
  gen$RANDFL <- as.factor(dplyr::case_when(
    gen$SAFFL == "Y" ~ "Y",
    .default = "N"
  ))
  gen$ITTFL <- as.factor(dplyr::case_when(
    gen$SAFFL == "Y" ~ "Y",
    .default = "N"
  ))
  gen$FASFL <- as.factor(dplyr::case_when(
    gen$SAFFL == "Y" ~ "Y",
    .default = "N"
  ))
  gen$PPROTFL <- as.factor(dplyr::case_when(
    gen$SAFFL == "Y" ~ "Y",
    .default = "N"
  ))
  gen$LSTSVDT <- dplyr::case_when(
    !is.na(gen$LSTALVDT) ~ gen$LSTALVDT,
    !is.na(gen$SCRFDT) ~ gen$SCRFDT
  )
  gen$EOSDY <- as.numeric(gen$EOSDT - gen$RANDDT + 1)
  gen$UNBLNDFL <- "Y"
  gen$RESCRNFL <- "Y"
  gen$DTHTRTFL <- dplyr::case_when(
    gen$TRTSDT <= gen$DTHDT & gen$DTHDT <= gen$TRTEDT + 30 ~ "Y",
    .default = NA
  )
  gen$DTHCAUSP <- dplyr::case_when(
    !is.na(gen$DTHDT) ~ "Death Cause Specify",
    .default = NA
  )
  gen$DTHAFTFL <- dplyr::case_when(
    gen$DTHDT > gen$TRTEDT ~ "Y",
    .default = NA
  )
  gen$DTH60TFL <- dplyr::case_when(
    gen$DTHDT <= gen$TRTSDT + 60 ~ "Y",
    .default = "N"
  )
  gen$UNBLNDDY <- as.numeric(dplyr::case_when(
    gen$UNBLNDFL == "Y" ~ gen$TRTSDT - gen$RANDDT + 1,
    .default = NA
  ))
  gen$UNBREAS <- dplyr::case_when(
    gen$UNBLNDFL == "Y" ~ "Unblind reason",
    .default = NA
  )
  gen$LDOSE <- as.numeric(20)
  gen$LDOSU <- "mg"
  gen$DTHTERM <- gen$DTHCAUS
  gen$LDSTODTH <- as.numeric(gen$DTHDT - gen$TRTEDT + 1)
  gen$DTHDY <- as.numeric(gen$DTHDT - gen$TRTSDT + 1)
  gen$DTHFL <- as.factor(gen$DTHFL)
  gen$DTH30FL <- as.factor(gen$DTH30FL)
  gen$AGEU <- as.factor(gen$AGEU)
  gen$EOSSTT <- as.factor(gen$EOSSTT)
  gen$REGION1 <- as.factor(gen$REGION1)
  gen$DTHA30FL <- as.factor(gen$DTHA30FL)
  gen$DTHB30FL <- as.factor(gen$DTHB30FL)

  # Add FASFL flag when TRT is NA
  gen$FASFL <- as.factor(dplyr::case_when(
    is.na(gen$TRT01P) ~ "N",
    .default = "Y"
  ))

  # remove NA TRTEDY
  gen <- dplyr::filter(gen, !is.na(TRTEDY))

  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    AGEGR1N = "Pooled Age Group 1 (N)",
    SEX_DECODE = "Sex",
    RACE_DECODE = "Race",
    ETHNIC_DECODE = "Ethnicity",
    WEIGHTBL = "Weight (kg)",
    WGTGR1N = "Weight Group 1 (N)",
    WGTGR1 = "Weight Group 1",
    HEIGHTBL = "Height (cm)",
    BSABL = "Body surface area (m2)",
    BMIBL = "Body mass index (kg/m2)",
    BMIBLG1N = "BMI at Baseline Group 1 (N)",
    BMIBLG1 = "BMI at Baseline Group 1",
    COUNTRY_DECODE = "Country",
    RFICDT = "Date of Informed Consent",
    RANDFL = "Randomized Flag",
    RACEGR1 = "Pooled Race Group 1",
    AGEGR1 = "Pooled Age Group 1",
    TRTEDY = "Treatment Relative End Day",
    RFICDT = "Date of Informed Consent",
    DTHTRTFL = "Death on Treatment Flag",
    DTHCAUSP = "Cause Spec for Death",
    LSTSVDT = "Last Subject Visit (SV) Date",
    EOTSTT = "End of Treatment Status",
    EOSDY = "Study Day of Study Termination",
    LDOSE = "Last Dose",
    LDOSU = "Last Dose Unit",
    AGEGR1N = "Pooled Age Group 1 (N)",
    RANUM = "Randomization Number",
    STRAT1R = "Strat Factor 1 Value Used for Rand",
    STRAT2R = "Strat Factor 2 Value Used for Rand",
    SCRNFL = "Screened Population Flag",
    DTHAFTFL = "Death After 30 Days of Last Treatment",
    DTH60TFL = "Death Within 60 Days of First Treatment",
    DTHTERM = "Reported Cause of Death",
    LDSTODTH = "Days from Last Dose to Death",
    RANDDTM = "Datetime of Randomization",
    FASFL = "Full Analysis Set Population Flag",
    ENRLFL = "Enrolled Population Flag",
    SCRFFL = "Screen Failure Flag",
    LSVISIT = "Last Study Visit",
    DCSREAS = "Reason for Discontinuation from Study",
    DCSREASP = "Reason Spec for Discont from Study",
    DCTREAS = "Reason for Discontinuation of Treatment",
    DCTREASP = "Reason Specify for Discont of Treatment",
    UNBLNDDY = "Study Day of Unblinding",
    UNBLNDFL = "Subject Blind Broken",
    UNBREAS = "Reason For Unblinding",
    DCSCREEN = "Reason for Discont During Screening",
    PPROTFL = "Per-Protocol Population Flag",
    LTVISIT = "Last Treatment Visit",
    DTHDY = "Study Day of Death",
    RESCRNFL = "Re-screened Flag",
    ITTFL = "Intent-To-Treat Population Flag"
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

adsl <- gen_adsl()
