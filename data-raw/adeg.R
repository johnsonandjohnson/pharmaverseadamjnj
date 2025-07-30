# Generate ADEG dataset

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(admiral)
library(formatters)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADEG dataset
gen_adeg <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get source data
  raw <- pharmaverseadam::adeg

  gen <- dplyr::filter(
    raw,
    (!is.na(PARAMCD)) &
      AVISIT != "Week 26" &
      (!is.na(AVISIT) & is.na(ATPT))
  )

  gen <- dplyr::mutate(
    gen,
    PARAMCD = as.factor(dplyr::case_when(
      PARAMCD == "HR" ~ "EGHRMN",
      PARAMCD == "QTLCR" ~ "PRAG",
      PARAMCD == "RRR" ~ "QRSAG",
      PARAMCD == "QT" ~ "QTC",
      PARAMCD == "QTCBR" ~ "QTCBAG",
      PARAMCD == "QTCFR" ~ "QTCFAG",
      PARAMCD == "RR" ~ "RRAG",
      PARAMCD == "EGINTP" ~ "INTP"
    )),
    PARAM = as.factor(dplyr::case_when(
      PARAMCD == "EGHRMN" ~ "ECG Mean Heart Rate (beats/min)",
      PARAMCD == "PRAG" ~ "PR Interval, Aggregate (msec)",
      PARAMCD == "QRSAG" ~ "QRS Duration, Aggregate (msec)",
      PARAMCD == "QTC" ~ "QT Interval, Corrected (msec)",
      PARAMCD == "QTCBAG" ~ "QTcB Interval, Aggregate (msec)",
      PARAMCD == "QTCFAG" ~ "QTcF Interval, Aggregate (msec)",
      PARAMCD == "RRAG" ~ "RR Interval, Aggregate (msec)",
      PARAMCD == "INTP" ~ "Interpretation"
    )),
    AVISIT = forcats::fct_reorder(
      as.factor(dplyr::case_when(
        AVISIT == "Baseline" ~ AVISIT,
        AVISIT == "Week 2" ~ "Month 1",
        AVISIT == "Week 4" ~ "Month 3",
        AVISIT == "Week 6" ~ "Month 6",
        AVISIT == "Week 8" ~ "Month 9",
        AVISIT == "Week 12" ~ "Month 12",
        AVISIT == "Week 16" ~ "Month 15",
        AVISIT == "Week 20" ~ "Month 18",
        AVISIT == "Week 24" ~ "Month 24"
      )),
      AVISITN
    ),
    ATPT = "BEFORE TREATMENT",
    ATPTN = 1,
    ADTM = as.POSIXct(paste0(ADT, "T06:00"), format = "%Y-%m-%dT%H:%M"),
    AVAL = dplyr::case_when(
      PARAMCD == "EGHRMN" ~
        as.numeric(sample(seq(0, 600), dplyr::n(), replace = TRUE)),
      PARAMCD == "PRAG" ~
        as.numeric(sample(seq(0, 600), dplyr::n(), replace = TRUE)),
      PARAMCD == "RRAG" ~
        as.numeric(sample(seq(0, 600), dplyr::n(), replace = TRUE)),
      PARAMCD == "QRSAG" ~
        as.numeric(sample(seq(0, 600), dplyr::n(), replace = TRUE)),
      PARAMCD == "QTC" ~
        as.numeric(sample(seq(200, 525), dplyr::n(), replace = TRUE)),
      PARAMCD == "QTCFAG" ~
        as.numeric(sample(seq(200, 525), dplyr::n(), replace = TRUE)),
      PARAMCD == "QTCBAG" ~
        as.numeric(sample(seq(200, 525), dplyr::n(), replace = TRUE)),
      .default = AVAL
    ),
    ABLFL = as.factor(dplyr::case_when(
      AVISIT == "Baseline" ~ "Y"
    )),
    TRTEMFL = dplyr::case_when(
      ONTRTFL == "Y" ~ "Y",
      .default = NA
    ),
    ANL01FL = "Y",
    ANL02FL = "Y",
    ANL03FL = dplyr::case_when(
      AVISIT == "Month 3" ~ "Y",
      .default = NA
    ),
    # ADY = as.numeric(ADT - TRTSDT + 1),
    AVALCAT1 = as.factor(dplyr::case_when(
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        AVAL <= 450 ~
        "<=450",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        AVAL > 450 &
        AVAL <= 480 ~
        ">450 to <=480 ",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        AVAL > 450 &
        AVAL <= 500 ~
        ">480 to <=500",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        AVAL > 500 ~
        ">500"
    )),
    APOBLFL = dplyr::case_when(
      AVISIT != "Baseline" ~ "Y",
      .default = NA
    ),
    CRIT1 = dplyr::case_when(
      PARAMCD == "EGHRMN" ~ "<50",
      PARAMCD == "PRAG" ~ "<120"
    ),
    CRIT1FL = dplyr::case_when(
      PARAMCD == "EGHRMN" & AVAL < 50 ~ "Y",
      PARAMCD == "PRAG" & AVAL < 120 ~ "Y"
    ),
    CRIT2 = dplyr::case_when(
      PARAMCD == "EGHRMN" ~ ">100",
      PARAMCD == "PRAG" ~ ">200"
    ),
    CRIT2FL = dplyr::case_when(
      PARAMCD == "EGHRMN" & AVAL > 100 ~ "Y",
      PARAMCD == "PRAG" & AVAL > 200 ~ "Y"
    ),
    EGCLSIG = as.factor(dplyr::case_when(
      AVAL > 500 & (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") ~ "Y",
      AVAL > 100 & PARAMCD == "EGHRMN" ~ "Y",
      AVAL > 200 & PARAMCD == "PRAG" ~ "Y",
      .default = "N"
    ))
  )

  # For INTP parameter, randomly populate AVALC with "NORMAL" and "ABNORMAL"
  # and set EGCLSIG accordingly
  intp_rows <- which(gen$PARAMCD == "INTP")
  # Randomly assign NORMAL or ABNORMAL to AVALC for INTP rows
  gen$AVALC[intp_rows] <- sample(c("NORMAL", "ABNORMAL"), length(intp_rows), replace = TRUE, prob = c(0.5, 0.5))

  # Set EGCLSIG based on AVALC
  normal_rows <- intp_rows[gen$AVALC[intp_rows] == "NORMAL"]
  abnormal_rows <- intp_rows[gen$AVALC[intp_rows] == "ABNORMAL"]

  # Set EGCLSIG to NA for NORMAL rows
  gen$EGCLSIG[normal_rows] <- NA

  # Randomly assign Y or N to EGCLSIG for ABNORMAL rows
  gen$EGCLSIG[abnormal_rows] <- sample(c("Y", "N"), length(abnormal_rows), replace = TRUE, prob = c(0.5, 0.5))

  # Convert EGCLSIG back to factor
  gen <- dplyr::mutate(
    gen,
    EGCLSIG = as.factor(EGCLSIG)
  )

  gen <- dplyr::select(gen, -BASE, -BNRIND)

  gen <- admiral::derive_var_base(
    gen,
    by_vars = exprs(USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  )

  gen <- admiral::derive_var_base(
    gen,
    by_vars = exprs(USUBJID, PARAMCD),
    source_var = ANRIND,
    new_var = BNRIND,
    filter = ABLFL == "Y"
  )

  gen <- admiral::derive_var_chg(gen)
  gen <- admiral::derive_var_pchg(gen)

  gen <- dplyr::mutate(
    gen,
    CHG = dplyr::case_when(
      ABLFL == "Y" ~ NA,
      .default = CHG
    ),
    PCHG = dplyr::case_when(
      ABLFL == "Y" ~ NA,
      .default = PCHG
    ),
    BASECAT1 = as.factor(dplyr::case_when(
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        BASE <= 450 ~
        "<=450",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        BASE > 450 &
        BASE <= 480 ~
        ">450 to <=480 ",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        BASE > 450 &
        BASE <= 500 ~
        ">480 to <=500",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        BASE > 500 ~
        ">500"
    )),
    CHGCAT1 = dplyr::case_when(
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        CHG <= 30 ~
        "<=30",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        CHG > 30 &
        CHG <= 60 ~
        ">30 to <=60",
      (PARAMCD == "QTC" | PARAMCD == "QTCBAG" | PARAMCD == "QTCFAG") &
        CHG > 60 ~
        ">60"
    )
  )

  source(file.path("data-raw", "adsl.R"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), "USUBJID")

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c(
    "TRT01A",
    "SAFFL",
    "STUDYID",
    "AGE",
    "SEX",
    "RACE_DECODE"
  )

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl %>%
    select(USUBJID, all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    message("Dropping shared vars from raw: ", paste(shared, collapse = ", "))
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID")

  # Define additional labels for variables added/modified or joined
  additional_labels <- list(
    PARAM = "Parameter",
    BASECAT1 = "Baseline Category 1",
    # Treatment flags and analysis flags
    TRTEMFL = "Treatment Emergent Analysis Flag",
    ANL02FL = "Analysis Flag 02-By Visit Value",
    ANL03FL = "Analysis Flag 03-Maximum Value",
    ANL01FL = "Analysis Flag 01-Analysis Value",
    APOBLFL = "Post-Baseline Record Flag",

    # Criteria fields
    CRIT1 = "Analysis Criterion 1",
    CRIT1FL = "Criterion 1 Evaluation Result Flag",
    CRIT2 = "Analysis Criterion 2",
    CRIT2FL = "Criterion 2 Evaluation Result Flag",
    EGCLSIG = "Clinically Significant",

    # Categorization variable
    BASECAT1 = "Baseline Category 1"
  )

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels using raw and additional_labels
  gen <- restore_labels(
    df = gen,
    orig_df = raw, # Use original pharmaverseadam::adeg
    additional_labels = additional_labels
  )

  return(gen)
}


adeg <- gen_adeg()
