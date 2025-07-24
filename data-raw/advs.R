# Generate ADVS dataset

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(formatters)
library(labelled)
library(admiral)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADVS dataset
gen_advs <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  raw <- pharmaverseadam::advs

  gen <- raw |>
    # nolint start
    dplyr::filter(
      (PARAMCD == "SYSBP" |
        PARAMCD == "DIABP" |
        PARAMCD == "PULSE" |
        PARAMCD == "TEMP" |
        PARAMCD == "WEIGHT") &
        DTYPE == "AVERAGE" &
        !is.na(AVISIT)
    ) |>
    # nolint end
    dplyr::mutate(
      AVALC = NA
    )

  # Create SYSBPO, DIABPO, and PULSEO Parameters
  gen_ortho <- gen |>
    dplyr::filter(
      PARAMCD == "SYSBP" | PARAMCD == "DIABP" | PARAMCD == "PULSE"
    ) |>
    dplyr::mutate(
      PARAMCD = dplyr::case_when(
        PARAMCD == "SYSBP" ~ "SYSBPO",
        PARAMCD == "DIABP" ~ "DIABPO",
        PARAMCD == "PULSE" ~ "PULSEO",
      ),
      PARAM = dplyr::case_when(
        PARAMCD == "SYSBPO" ~ "Systolic Blood Pressure Orthostatic (mmHg)",
        PARAMCD == "DIABPO" ~ "Diastolic Blood Pressure Orthostatic (mmHg)",
        PARAMCD == "PULSEO" ~ "Pulse Rate Orthostatic (beats/min)",
      ),
      AVAL = dplyr::case_when(
        PARAMCD == "SYSBPO" ~
          as.numeric(sample(seq(-30, 30), dplyr::n(), replace = TRUE)),
        PARAMCD == "DIABPO" ~
          as.numeric(sample(seq(-25, 25), dplyr::n(), replace = TRUE)),
        PARAMCD == "PULSEO" ~
          as.numeric(sample(seq(-25, 25), dplyr::n(), replace = TRUE)),
      )
    )

  gen_ortho_der <- gen_ortho |>
    dplyr::filter(PARAMCD == "SYSBPO" | PARAMCD == "DIABPO") |>
    dplyr::mutate(
      PARAMCD = dplyr::case_when(
        PARAMCD == "SYSBPO" ~ "ORTHYPS",
        PARAMCD == "DIABPO" ~ "ORTHYPD",
      ),
      PARAM = dplyr::case_when(
        PARAMCD == "ORTHYPS" ~ "SBP (STD-SUP) <-20",
        PARAMCD == "ORTHYPD" ~ "DBP (STD-SUP) <-10",
      ),
      AVALC = dplyr::case_when(
        PARAMCD == "ORTHYPS" & AVAL < -20 ~ "Y",
        PARAMCD == "ORTHYPS" ~ "N",
        PARAMCD == "ORTHYPD" & AVAL < -10 ~ "Y",
        PARAMCD == "ORTHYPD" ~ "N",
      ),
      AVAL = NA,
    )

  gen_orthyps <- gen_ortho_der |>
    dplyr::filter(PARAMCD == "ORTHYPS") |>
    dplyr::mutate(
      ORTHYPS = AVALC
    ) |>
    dplyr::select(-AVALC)

  gen_orthypd <- gen_ortho_der |>
    dplyr::filter(PARAMCD == "ORTHYPD") |>
    dplyr::mutate(
      ORTHYPD = AVALC
    ) |>
    dplyr::select(STUDYID, USUBJID, AVISITN, ORTHYPD)

  gen_orthyp <- dplyr::inner_join(
    gen_orthyps,
    gen_orthypd,
    by = c("STUDYID", "USUBJID", "AVISITN"),
    copy = FALSE,
    suffix = c(".x", ".y"),
    keep = FALSE,
    na_matches = "na"
  ) |>
    dplyr::mutate(
      PARAMCD = "ORTHYP",
      PARAM = "Orthostatic Hypotension",
      AVALC = dplyr::case_when(
        ORTHYPS == "Y" |
          ORTHYPD == "Y" ~
          "Y",
        .default = "N"
      )
    ) |>
    dplyr::select(-ORTHYPS, -ORTHYPD)

  gen <- rbind(gen, gen_ortho, gen_ortho_der, gen_orthyp) |>
    dplyr::mutate(
      ABLFL = dplyr::case_when(
        AVISIT == "Baseline" ~ "Y",
        .default = NA
      ),
      ATPT = "BEFORE TREATMENT",
      ATPTN = 1,
      ANL01FL = "Y",
      ANL02FL = "Y",
      APOBLFL = dplyr::case_when(
        AVISIT == "Baseline" ~ NA,
        .default = "Y"
      ),
      ONTRTFL = dplyr::case_when(
        AVISIT == "Baseline" ~ NA,
        .default = "Y"
      ),
      AVALCAT1 = dplyr::case_when(
        PARAMCD == "SYSBP" & AVAL < 90 ~ "<90",
        PARAMCD == "SYSBP" & AVAL >= 90 & AVAL < 119 ~ ">=90 to 119",
        PARAMCD == "SYSBP" & AVAL >= 120 & AVAL < 139 ~ ">=120 to 139",
        PARAMCD == "SYSBP" & AVAL >= 140 & AVAL < 159 ~ ">=140 to 159",
        PARAMCD == "SYSBP" & AVAL >= 160 & AVAL < 179 ~ ">=160 to 179",
        PARAMCD == "SYSBP" & AVAL >= 180 ~ ">=180",
        PARAMCD == "DIABP" & AVAL < 60 ~ "<60",
        PARAMCD == "DIABP" & AVAL >= 60 & AVAL < 89 ~ ">=60 to 89",
        PARAMCD == "DIABP" & AVAL >= 90 & AVAL < 109 ~ ">=90 to 109",
        PARAMCD == "DIABP" & AVAL >= 110 & AVAL < 119 ~ ">=110 to 119",
        PARAMCD == "DIABP" & AVAL >= 120 ~ ">=120",
      ),
      AVALCA1N = dplyr::case_when(
        PARAMCD == "SYSBP" & AVAL < 90 ~ 1,
        PARAMCD == "SYSBP" & AVAL >= 90 & AVAL < 119 ~ 2,
        PARAMCD == "SYSBP" & AVAL >= 120 & AVAL < 139 ~ 3,
        PARAMCD == "SYSBP" & AVAL >= 140 & AVAL < 159 ~ 4,
        PARAMCD == "SYSBP" & AVAL >= 160 & AVAL < 179 ~ 5,
        PARAMCD == "SYSBP" & AVAL >= 180 ~ 6,
        PARAMCD == "DIABP" & AVAL < 60 ~ 1,
        PARAMCD == "DIABP" & AVAL >= 60 & AVAL < 89 ~ 2,
        PARAMCD == "DIABP" & AVAL >= 90 & AVAL < 109 ~ 3,
        PARAMCD == "DIABP" & AVAL >= 110 & AVAL < 119 ~ 4,
        PARAMCD == "DIABP" & AVAL >= 120 ~ 5,
      ),
    ) |>
    dplyr::select(-BASE, -BNRIND) |>
    admiral::derive_var_base(
      by_vars = exprs(USUBJID, PARAMCD),
      source_var = AVAL,
      new_var = BASE,
      filter = ABLFL == "Y"
    ) |>
    admiral::derive_var_base(
      by_vars = exprs(USUBJID, PARAMCD),
      source_var = ANRIND,
      new_var = BNRIND,
      filter = ABLFL == "Y"
    ) |>
    admiral::derive_var_chg() |>
    admiral::derive_var_pchg() |>
    dplyr::mutate(
      CHG = dplyr::case_when(
        ABLFL == "Y" ~ NA,
        .default = CHG
      ),
      PCHG = dplyr::case_when(
        ABLFL == "Y" ~ NA,
        .default = PCHG
      ),
      ADTM = as.POSIXct(paste0(ADT, "T06:00"), format = "%Y-%m-%dT%H:%M"),
      CRIT1 = dplyr::case_when(
        PARAMCD == "SYSBP" ~
          "<90 mmHg and with >30 mmHg decrease from baseline",
        PARAMCD == "DIABP" ~
          "<50 mmHg and with >20 mmHg decrease from baseline",
        PARAMCD == "PULSE" ~ "<50 bpm and with >20 bpm decrease from baseline",
        PARAMCD == "WEIGHT" ~ "decrease 10% kg from baseline",
        PARAMCD == "RESP" ~ ">20 breaths per minute",
        PARAMCD == "TEMP" ~ ">38 and with >=1 increase from baseline",
        PARAMCD == "PULSEO" ~ "Orthostatic Pulse Rate >15",
        PARAMCD == "SYSBPO" ~ "Orthostatic SBP <-20",
        PARAMCD == "DIABPO" ~ "Orthostatic DBP <-10",
      ),
      CRIT1FL = dplyr::case_when(
        PARAMCD == "SYSBP" & AVAL < 90 & CHG < -30 ~ "Y",
        PARAMCD == "SYSBP" ~ "N",
        PARAMCD == "DIABP" & AVAL < 50 & CHG < -20 ~ "Y",
        PARAMCD == "DIABP" ~ "N",
        PARAMCD == "PULSE" & AVAL < 50 & CHG < -20 ~ "Y",
        PARAMCD == "PULSE" ~ "N",
        PARAMCD == "WEIGHT" & PCHG < -10 ~ "Y",
        PARAMCD == "WEIGHT" ~ "N",
        PARAMCD == "RESP" & AVAL > 20 ~ "Y",
        PARAMCD == "RESP" ~ "N",
        PARAMCD == "TEMP" & AVAL > 38 & CHG >= 1 ~ "Y",
        PARAMCD == "TEMP" ~ "N",
        PARAMCD == "PULSEO" & AVAL > 15 ~ "Y",
        PARAMCD == "PULSEO" ~ "N",
        PARAMCD == "SYSBPO" & AVAL < -20 ~ "Y",
        PARAMCD == "SYSBPO" ~ "N",
        PARAMCD == "DIABPO" & AVAL < -10 ~ "Y",
        PARAMCD == "DIABPO" ~ "N",
      ),
      CRIT2 = dplyr::case_when(
        PARAMCD == "SYSBP" ~
          ">180 mmHg and with >40 mmHg increase from baseline",
        PARAMCD == "DIABP" ~
          ">105 mmHg and with >30 mmHg increase from baseline",
        PARAMCD == "PULSE" ~ ">120 bpm and with >30 bpm increase from baseline",
        PARAMCD == "WEIGHT" ~ "increase 10% kg from baseline",
      ),
      CRIT2FL = dplyr::case_when(
        PARAMCD == "SYSBP" & AVAL > 180 & CHG > 40 ~ "Y",
        PARAMCD == "SYSBP" ~ "N",
        PARAMCD == "DIABP" & AVAL > 105 & CHG > 30 ~ "Y",
        PARAMCD == "DIABP" ~ "N",
        PARAMCD == "PULSE" & AVAL > 120 & CHG > 30 ~ "Y",
        PARAMCD == "PULSE" ~ "N",
        PARAMCD == "WEIGHT" & CHG > 10 ~ "Y",
        PARAMCD == "WEIGHT" ~ "N",
      ),
      CRIT3 = dplyr::case_when(
        PARAMCD == "SYSBP" ~ "Systolic blood pressure<90",
        PARAMCD == "DIABP" ~ "Diastolic blood pressure<60",
      ),
      CRIT3FL = dplyr::case_when(
        PARAMCD == "SYSBP" & AVAL < 90 ~ "Y",
        PARAMCD == "SYSBP" ~ "N",
        PARAMCD == "DIABP" & AVAL < 60 ~ "Y",
        PARAMCD == "DIABP" ~ "N",
      ),
      ATOXDSCL = dplyr::case_when(
        PARAMCD == "SYSBP" ~ "Hypotension (systolic)",
        PARAMCD == "PULSE" ~ "Bradycardia",
      ),
      ATOXDSCH = dplyr::case_when(
        PARAMCD == "SYSBP" ~ "Hypertension (systolic)",
        PARAMCD == "PULSE" ~ "Tachycardia",
        PARAMCD == "DIABP" ~ "Hypertension (diastolic)",
        PARAMCD == "RESP" ~ "Respiratory Rate",
        PARAMCD == "TEMP" ~ "Fever",
      ),
      ATOXGRL = dplyr::case_when(
        PARAMCD == "PULSE" & AVAL >= 55 & AVAL < 999999 ~ 0,
        PARAMCD == "PULSE" & AVAL >= 50 & AVAL < 55 ~ 1,
        PARAMCD == "PULSE" & AVAL >= 45 & AVAL < 50 ~ 2,
        PARAMCD == "PULSE" & AVAL >= -999999 & AVAL < 45 ~ 3,
        PARAMCD == "SYSBP" & AVAL >= 90 & AVAL < 999999 ~ 0,
        PARAMCD == "SYSBP" & AVAL >= 85 & AVAL < 90 ~ 1,
        PARAMCD == "SYSBP" & AVAL >= 80 & AVAL < 85 ~ 2,
        PARAMCD == "SYSBP" & AVAL >= -999999 & AVAL < 80 ~ 3,
      ),
      ATOXGRH = dplyr::case_when(
        PARAMCD == "SYSBP" & AVAL > -999999 & AVAL <= 140 ~ 0,
        PARAMCD == "SYSBP" & AVAL > 140 & AVAL <= 150 ~ 1,
        PARAMCD == "SYSBP" & AVAL > 150 & AVAL <= 155 ~ 2,
        PARAMCD == "SYSBP" & AVAL > 155 & AVAL < 99999999 ~ 3,
        PARAMCD == "PULSE" & AVAL > -999999 & AVAL <= 100 ~ 0,
        PARAMCD == "PULSE" & AVAL > 100 & AVAL <= 115 ~ 1,
        PARAMCD == "PULSE" & AVAL > 115 & AVAL <= 130 ~ 2,
        PARAMCD == "PULSE" & AVAL > 130 & AVAL < 99999999 ~ 3,
        PARAMCD == "DIABP" & AVAL > -999999 & AVAL <= 90 ~ 0,
        PARAMCD == "DIABP" & AVAL > 90 & AVAL <= 95 ~ 1,
        PARAMCD == "DIABP" & AVAL > 95 & AVAL <= 100 ~ 2,
        PARAMCD == "DIABP" & AVAL > 100 & AVAL < 99999999 ~ 3,
        PARAMCD == "RESP" & AVAL > -999999 & AVAL < 17 ~ 0,
        PARAMCD == "RESP" & AVAL >= 17 & AVAL <= 20 ~ 1,
        PARAMCD == "RESP" & AVAL > 20 & AVAL <= 25 ~ 2,
        PARAMCD == "RESP" & AVAL > 25 & AVAL < 99999999 ~ 3,
        PARAMCD == "TEMP" & AVAL > -999999 & AVAL < 38.0 ~ 0,
        PARAMCD == "TEMP" & AVAL >= 38.0 & AVAL < 38.5 ~ 1,
        PARAMCD == "TEMP" & AVAL >= 38.5 & AVAL < 39.0 ~ 2,
        PARAMCD == "TEMP" & AVAL >= 39.0 & AVAL <= 40.0 ~ 3,
        PARAMCD == "TEMP" & AVAL > 40.0 & AVAL < 99999999 ~ 4,
      ),
      AVISITN = dplyr::case_when(
        AVISIT == "Baseline" ~ 1,
        AVISIT == "Week 2" ~ 2,
        AVISIT == "Week 4" ~ 3,
        AVISIT == "Week 6" ~ 4,
        AVISIT == "Week 8" ~ 5,
        AVISIT == "Week 12" ~ 6,
        AVISIT == "Week 16" ~ 7,
        AVISIT == "Week 20" ~ 8,
        AVISIT == "Week 24" ~ 9,
        AVISIT == "Week 26" ~ 23,
      ),
      AVISIT = forcats::fct_reorder(
        as.factor(dplyr::case_when(
          AVISIT == "Baseline" ~ "Screening",
          AVISIT == "Week 2" ~ "Cycle 02",
          AVISIT == "Week 4" ~ "Cycle 03",
          AVISIT == "Week 6" ~ "Cycle 04",
          AVISIT == "Week 8" ~ "Cycle 05",
          AVISIT == "Week 12" ~ "Cycle 06",
          AVISIT == "Week 16" ~ "Cycle 07",
          AVISIT == "Week 20" ~ "Cycle 08",
          AVISIT == "Week 24" ~ "Cycle 09",
          AVISIT == "Week 26" ~ "End Of Treatment",
        )),
        AVISITN
      ),
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ATOXGR = factor(pmax(ATOXGRL, ATOXGRH, na.rm = TRUE), levels = c(0, 1, 2, 3, 4)),
      ATOXGRL = as.factor(ATOXGRL),
      ATOXGRH = as.factor(ATOXGRH)
    ) |> # recreate factors
    dplyr::ungroup()

  gen_add_avisit1 <- gen |>
    dplyr::filter(AVISIT != "Screening") |>
    dplyr::mutate(
      AVISITN = dplyr::case_when(
        AVISIT == "Cycle 02" ~ 10,
        AVISIT == "Cycle 03" ~ 11,
        AVISIT == "Cycle 04" ~ 12,
        AVISIT == "Cycle 05" ~ 13,
        AVISIT == "Cycle 06" ~ 14,
        AVISIT == "Cycle 07" ~ 15,
        AVISIT == "Cycle 08" ~ 16,
        AVISIT == "Cycle 09" ~ 17,
        AVISIT == "End Of Treatment" ~ 18,
      ),
      AVISIT = dplyr::case_when(
        AVISIT == "Cycle 02" ~ "Cycle 10",
        AVISIT == "Cycle 03" ~ "Cycle 11",
        AVISIT == "Cycle 04" ~ "Cycle 12",
        AVISIT == "Cycle 05" ~ "Cycle 13",
        AVISIT == "Cycle 06" ~ "Cycle 15",
        AVISIT == "Cycle 07" ~ "Cycle 17",
        AVISIT == "Cycle 08" ~ "Cycle 19",
        AVISIT == "Cycle 09" ~ "Cycle 21",
        AVISIT == "End Of Treatment" ~ "Cycle 23",
      ),
      ADT = ADT + 365,
    )

  gen_add_avisit2 <- gen |>
    dplyr::filter(
      AVISIT == "Cycle 08" | AVISIT == "Cycle 09" | AVISIT == "End Of Treatment"
    ) |>
    dplyr::mutate(
      AVISITN = dplyr::case_when(
        AVISIT == "Cycle 08" ~ 19,
        AVISIT == "Cycle 09" ~ 21,
        AVISIT == "End Of Treatment" ~ 22
      ),
      AVISIT = dplyr::case_when(
        AVISIT == "Cycle 08" ~ "Cycle 22",
        AVISIT == "Cycle 09" ~ "Cycle 25",
        AVISIT == "End Of Treatment" ~ "Cycle 29"
      ),
      ADT = ADT + 730,
    )

  gen <- rbind(gen, gen_add_avisit1, gen_add_avisit2) |>
    dplyr::mutate(
      AVISIT = forcats::fct_reorder(as.factor(AVISIT), AVISITN),
    ) |>
    admiral::restrict_derivation(
      derivation = admiral::derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(desc(AVAL), ADT),
        new_var = ANL06FL,
        true_value = "Y",
        false_value = NA,
        mode = "last"
      ),
      filter = AVISIT != "Screening"
    ) |>
    admiral::restrict_derivation(
      derivation = admiral::derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(ATOXGRH, ADT),
        new_var = ANL05FL,
        true_value = "Y",
        false_value = NA,
        mode = "last"
      ),
      filter = AVISIT != "Screening"
    ) |>
    admiral::restrict_derivation(
      derivation = admiral::derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(desc(ATOXGRL), ADT),
        new_var = ANL04FL,
        true_value = "Y",
        false_value = NA,
        mode = "last"
      ),
      filter = AVISIT != "Screening"
    ) |>
    admiral::restrict_derivation(
      derivation = admiral::derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, PARAMCD),
        order = exprs(AVAL, ADT),
        new_var = ANL03FL,
        true_value = "Y",
        false_value = NA,
        mode = "last"
      ),
      filter = AVISIT != "Screening"
    ) |>
    mutate(
      TRTEMFL = dplyr::case_when(
        ADT > TRTSDT & ADT < TRTEDT ~ "Y",
        .default = "N"
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

  additional_labels <- list(
    SAFFL = "Safety Population Flag",
    AVALC = "Analysis Value (C)",
    AVALU = "Analysis Value Unit",
    ANL01FL = "Analysis Flag 01",
    COUNTRY_DECODE = "Country",
    PARCAT1 = "Parameter Category 1",
    PARCAT2 = "Parameter Category 2",
    ADTM = "Analysis Date/Time",
    ATPT = "Analysis Time Point",
    TR01SDT = "Treatment Start Date",
    TR01EDT = "Treatment End Date",
    ANL02FL = "Analysis Flag 02-By Visit Value",
    APOBLFL = "Post-Baseline Record Flag",
    BASE = "Baseline Value",
    BNRIND = "Baseline Reference Range Indicator",
    CRIT1 = "Analysis Criterion 1",
    CRIT1FL = "Criterion 1 Evaluation Result Flag",
    CRIT2 = "Analysis Criterion 2",
    CRIT2FL = "Criterion 2 Evaluation Result Flag",
    CRIT3 = "Analysis Criterion 3",
    CRIT3FL = "Criterion 3 Evaluation Result Flag",
    ATOXDSCL = "Analysis Toxicity Description Low",
    ATOXDSCH = "Analysis Toxicity Description High",
    ATOXGRL = "Analysis Toxicity Grade Low",
    ATOXGRH = "Analysis Toxicity Grade High",
    ATOXGR = "Analysis Toxicity Grade",
    ANL06FL = "Analysis Flag 06-Minimum Value",
    ANL05FL = "Analysis Flag 05-Worst Tox Grade High",
    ANL04FL = "Analysis Flag 04-Worst Value",
    ANL03FL = "Analysis Flag 03-Maximum Value",
    TRTEMFL = "Treatment Emergent Analysis Flag"
  )

  # Restore labels and handle NA
  gen <- df_na(gen, char_as_factor = TRUE)

  gen <- restore_labels(
    df = gen,
    orig_df = raw,
    additional_labels = additional_labels
  )

  return(gen)
}

# Generate and save the dataset
advs <- gen_advs()
