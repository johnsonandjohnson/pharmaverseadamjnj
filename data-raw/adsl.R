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

  gen <- dplyr::mutate(
    raw,
    RACE = factor(
      sample(
        c(
          "AMERICAN INDIAN OR ALASKA NATIVE",
          "ASIAN",
          "BLACK OR AFRICAN AMERICAN",
          "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
          "WHITE",
          "MULTIPLE",
          "NOT REPORTED",
          "UNKNOWN",
          "OTHER"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      levels = c(
        "AMERICAN INDIAN OR ALASKA NATIVE",
        "ASIAN",
        "BLACK OR AFRICAN AMERICAN",
        "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
        "WHITE",
        "MULTIPLE",
        "NOT REPORTED",
        "UNKNOWN",
        "OTHER"
      )
    ),
    ETHNIC = factor(
      sample(
        c(
          "HISPANIC OR LATINO",
          "NOT HISPANIC OR LATINO",
          "UNKNOWN",
          "NOT REPORTED"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      levels = c(
        "HISPANIC OR LATINO",
        "NOT HISPANIC OR LATINO",
        "UNKNOWN",
        "NOT REPORTED"
      )
    )
  )
  gen$TRT01P <- as.factor(gen$TRT01P)
  gen$TRT01P <- droplevels(as.factor(dplyr::case_when(
    gen$TRT01P == "Screen Failure" ~ NA,
    .default = gen$TRT01P
  )))
  gen$TRT01PN <- dplyr::case_when(
    gen$TRT01P == "Xanomeline Low Dose" ~ 1,
    gen$TRT01P == "Xanomeline High Dose" ~ 2,
    gen$TRT01P == "Placebo" ~ 3
  )
  gen$TRT01P <- forcats::fct_reorder(gen$TRT01P, gen$TRT01PN, .na_rm = TRUE)
  gen$TRT01A <- as.factor(gen$TRT01A)
  gen$TRT01A <- droplevels(as.factor(dplyr::case_when(
    gen$TRT01A == "Screen Failure" ~ NA,
    .default = gen$TRT01A
  )))
  gen$TRT01AN <- dplyr::case_when(
    gen$TRT01A == "Xanomeline Low Dose" ~ 1,
    gen$TRT01A == "Xanomeline High Dose" ~ 2,
    gen$TRT01A == "Placebo" ~ 3
  )
  gen$TRT01A <- forcats::fct_reorder(gen$TRT01A, gen$TRT01AN, .na_rm = TRUE)

  gen$AGEGR1 <- factor(
    dplyr::case_when(
      gen$AGE >= 18 & gen$AGE < 65 ~ ">=18 to <65",
      gen$AGE >= 65 & gen$AGE < 75 ~ ">=65 to <75",
      gen$AGE >= 75 ~ ">=75"
    ),
    levels = c(
      ">=18 to <65",
      ">=65 to <75",
      ">=75"
    )
  )
  gen$AGEGR1N <- dplyr::case_when(
    gen$AGEGR1 == ">=18 to <65" ~ 1,
    gen$AGEGR1 == ">=65 to <75" ~ 2,
    gen$AGEGR1 == ">=75" ~ 3
  )

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
      gen$WEIGHTBL >= 30 & gen$WEIGHTBL < 60 ~ ">=30 to <60",
      gen$WEIGHTBL >= 60 & gen$WEIGHTBL < 90 ~ ">=60 to <90",
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
      gen$BMIBL >= 18.5 & gen$BMIBL < 25 ~ "Normal >=18.5 to <25",
      gen$BMIBL >= 25 & gen$BMIBL < 30 ~ "Overweight >=25 to <30",
      gen$BMIBL >= 30 ~ "Obese >=30"
    )),
    gen$BMIBLG1N
  )
  gen$SEX <- as.factor(gen$SEX)
  gen$REGION1 <- "North America"
  gen$RACEGR1 <- as.factor(gen$RACEGR1)
  gen$RFICDTC <- gen$DMDTC
  gen$RFICDT <- as.Date(gen$DMDTC)


  # Stratification factors
  gen$STRAT1D <- as.factor("Description of Stratification Factor 1")
  gen$STRAT2D <- as.factor("Description of Stratification Factor 2")

  # Randomized values (exactly ~50% each, shuffled for randomness)
  n_subj <- nrow(gen)
  n_first <- floor(n_subj / 2)
  vals1 <- c(
    rep("Stratification Factor 1 Value 1", n_first),
    rep("Stratification Factor 1 Value 2", n_subj - n_first)
  )
  vals2 <- c(
    rep("Stratification Factor 2 Value 1", n_first),
    rep("Stratification Factor 2 Value 2", n_subj - n_first)
  )
  gen$STRAT1R <- factor(sample(vals1),
    levels = c(
      "Stratification Factor 1 Value 1",
      "Stratification Factor 1 Value 2"
    )
  )
  gen$STRAT2R <- factor(sample(vals2),
    levels = c(
      "Stratification Factor 2 Value 1",
      "Stratification Factor 2 Value 2"
    )
  )
  gen$RANUM <- as.factor("1000001")
  gen$RANDDTM <- as.POSIXct(
    paste0(gen$RANDDT, " 11:59"),
    format = "%Y-%m-%d %H:%M"
  )
  # Randomly assign "CONTINUING" to some subjects for EOTSTT & EOSSTT
  set.seed(seed + 1)

  gen$EOTSTT <- factor(
    dplyr::case_when(
      sample(c(TRUE, FALSE), nrow(gen), replace = TRUE, prob = c(0.2, 0.8)) ~ "ONGOING",
      .default = as.character(gen$EOSSTT)
    ),
    levels = c("ONGOING", "COMPLETED", "DISCONTINUED")
  )

  gen$EOSSTT <- factor(
    dplyr::case_when(
      sample(c(TRUE, FALSE), nrow(gen), replace = TRUE, prob = c(0.2, 0.8)) ~ "ONGOING",
      .default = as.character(gen$EOSSTT)
    ),
    levels = c("ONGOING", "COMPLETED", "DISCONTINUED")
  )

  gen$DCTREAS <- factor(
    dplyr::case_when(
      gen$EOTSTT == "DISCONTINUED" ~ "Other",
      .default = NA
    ),
    levels = c("Other")
  )
  gen$LTVISIT <- as.factor("Last Treatment Visit")
  gen$DCTREASP <- dplyr::case_when(
    gen$DCTREAS == "Other" ~ "specify text",
    .default = NA
  )
  gen$DCTDT <- dplyr::case_when(
    gen$EOTSTT == "DISCONTINUED" ~ gen$EOSDT,
    .default = NA
  )

  gen$DCSREAS <- factor(
    dplyr::case_when(
      gen$EOSSTT == "DISCONTINUED" ~ "Other",
      .default = NA
    ),
    levels = c("Other")
  )
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

  gen$SCRFFL <- factor(
    dplyr::case_when(
      !is.na(gen$SCRFDT) ~ "Y",
      .default = "N"
    ),
    levels = c("Y", "N")
  )

  gen$DCSCREEN <- factor(
    dplyr::case_when(
      !is.na(gen$SCRFDT) ~ "Failure to meet eligibility criteria",
      .default = NA
    ),
    levels = c("Failure to meet eligibility criteria")
  )
  gen$SAFFL <- dplyr::case_when(
    is.na(gen$SAFFL) ~ "N",
    # Randomly set SAFFL="N" for 10% of subjects with non-NA TRT01P
    !is.na(gen$TRT01P) & sample(c(TRUE, FALSE), nrow(gen), replace = TRUE, prob = c(0.3, 0.7)) ~ "N",
    .default = gen$SAFFL
  )

  gen$ENRLFL <- factor(
    dplyr::case_when(
      gen$SAFFL == "Y" ~ "Y",
      .default = "N"
    ),
    levels = c("Y", "N")
  )
  gen$RANDFL <- ifelse(gen$SAFFL == "Y", "Y", "N")
  gen$RANDFL <- as.factor(gen$RANDFL)

  gen$ITTFL <- factor(
    dplyr::case_when(
      gen$SAFFL == "Y" ~ "Y",
      .default = "N"
    ),
    levels = c("Y", "N")
  )

  gen$FASFL <- factor(
    dplyr::case_when(
      gen$SAFFL == "Y" ~ "Y",
      .default = "N"
    ),
    levels = c("Y", "N")
  )

  gen$PPROTFL <- factor(
    dplyr::case_when(
      gen$SAFFL == "Y" ~ "Y",
      .default = "N"
    ),
    levels = c("Y", "N")
  )
  gen$LSTSVDT <- dplyr::case_when(
    !is.na(gen$LSTALVDT) ~ gen$LSTALVDT,
    !is.na(gen$SCRFDT) ~ gen$SCRFDT
  )
  gen$LASTCTDT <- dplyr::case_when(
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
  gen$DTHB60FL <- dplyr::case_when(
    gen$DTHDT <= gen$TRTSDT + 60 ~ "Y",
    .default = "N"
  )
  gen$UNBLNDDT <- as.Date(dplyr::case_when(
    gen$UNBLNDFL == "Y" ~ gen$TRTSDT + 1,
    .default = NA
  ))
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

  gen$FASFL <- factor(
    dplyr::case_when(
      is.na(gen$TRT01P) ~ "N",
      .default = "Y"
    ),
    levels = c("Y", "N")
  )

  # remove NA TRTEDY
  gen <- dplyr::filter(gen, !is.na(TRTEDY))

  gen <- dplyr::mutate(
    gen,
    PKFL = dplyr::if_else(USUBJID %in% pharmaverseadam::adpc$USUBJID, "Y", "N"),
    DIABETFL = sample(c("N", "Y"), NROW(gen), TRUE, c(0.8, 0.2))
  )

  gen <- dplyr::mutate(
    gen,
    IMFL = PKFL
  )

  gen <- dplyr::mutate(
    gen,
    DCTADY = as.numeric(DCTDT - TRTSDT + 1)
  )

  gen <- dplyr::mutate(
    gen,
    SAFEXRS = dplyr::case_when(
      toupper(SAFFL) != "Y" ~ "Exclusion reason for safety analysis set",
      .default = NA
    ),
    FASEXRS = dplyr::case_when(
      toupper(FASFL) != "Y" ~ "Exclusion reason for full analysis set",
      .default = NA
    ),
    PPREXRS = dplyr::case_when(
      toupper(PPROTFL) != "Y" ~ "Exclusion reason for per-protocol analysis set",
      .default = NA
    ),
    PKEXRES = dplyr::case_when(
      toupper(PKFL) != "Y" ~ "Exclusion reason for pharmacokinetics analysis set",
      .default = NA
    ),
    IMEXRES = dplyr::case_when(
      toupper(IMFL) != "Y" ~ "Exclusion reason for immunogenicity analysis set",
      .default = NA
    ),
  )

  gen <- dplyr::mutate(
    gen,
    DCSCREEN = case_when(
      USUBJID == "01-701-1240" ~ "Subject refused to sign informed consent",
      .default = DCSCREEN
    ),
    RESCRNFL = if_else(
      SCRFFL == "Y" & runif(n()) < 0.5, "Y", NA_character_
    )
  )


  gen <- gen |>
    dplyr::mutate(
      COHORT = factor(
        dplyr::case_match(
          ARM,
          "Placebo" ~ "Cohort 1",
          "Xanomeline High Dose" ~ "Cohort 2",
          "Xanomeline Low Dose" ~ "Cohort 3",
          "Screen Failure" ~ NA_character_,
          .default = NA_character_
        ),
        levels = c("Cohort 1", "Cohort 2", "Cohort 3")
      ),
      GROUP = factor(
        dplyr::case_match(
          ARM,
          "Placebo" ~ "Group 1",
          "Xanomeline High Dose" ~ "Group 2",
          "Xanomeline Low Dose" ~ "Group 3",
          "Screen Failure" ~ NA_character_,
          .default = NA_character_
        ),
        levels = c("Group 1", "Group 2", "Group 3")
      ),
      EOTDT = TRTEDT
    )

  # Define additional labels for new variables not in source dataset
  additional_labels <- list(
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    AGEGR1N = "Pooled Age Group 1 (N)",
    WEIGHTBL = "Weight (kg)",
    WGTGR1N = "Weight Group 1 (N)",
    WGTGR1 = "Weight Group 1",
    HEIGHTBL = "Height (cm)",
    BSABL = "Body surface area (m2)",
    BMIBL = "Body mass index (kg/m2)",
    BMIBLG1N = "BMI at Baseline Group 1 (N)",
    BMIBLG1 = "BMI at Baseline Group 1",
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
    STRAT1D = "Description of Stratification Factor 1",
    STRAT1R = "Strat Factor 1 Value Used for Rand",
    STRAT2D = "Description of Stratification Factor 2",
    STRAT2R = "Strat Factor 2 Value Used for Rand",
    SCRNFL = "Screened Population Flag",
    DTHAFTFL = "Death After 30 Days of Last Treatment",
    DTHB60FL = "Death Within 60 Days of First Treatment",
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
    UNBLNDDT = "Date of Unblinding",
    UNBLNDFL = "Subject Blind Broken",
    UNBREAS = "Reason For Unblinding",
    DCSCREEN = "Reason for Discont During Screening",
    PPROTFL = "Per-Protocol Population Flag",
    LTVISIT = "Last Treatment Visit",
    LASTCTDT = "Last Contact Date",
    DTHDY = "Study Day of Death",
    RESCRNFL = "Re-screened Flag",
    ITTFL = "Intent-To-Treat Population Flag",
    PKFL = "Pharmacokinetic Population Flag",
    IMFL = "Immunogenicity Population Flag",
    DIABETFL = "History of Diabetes",
    DCTADY = "Study Day of Treatment Discontinuation",
    SAFEXRS = "Reason for Excl from Safety Population",
    FASEXRS = "Reason for Excl from Full Analysis Set",
    PPREXRS = "Reason for Excl from Per-Prot Population",
    PKEXRES = "Reason for Excl from Pharmacokin Pop",
    IMEXRES = "Reason for Excl from Immunogen Pop",
    COHORT = "Cohort",
    GROUP = "Analysis Group",
    EOTDT = "End-of-Treatment Date",
    BRTHDTC = "Date/Time of Birth"
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
