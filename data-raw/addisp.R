# Generate ADDISP dataset

# Load necessary libraries
library(dplyr)
library(tidyr)
library(pharmaverseadam)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADDISP dataset based on addisp.md
gen_addisp <- function(seed = 123) {
  set.seed(seed)

  # Source related datasets
  source(file.path("data-raw", "adsl.R"))

  adsl <- adsl |>
    mutate(
      DCTDT = if (!("DCTDT" %in% names(pick(everything())))) as.Date(NA) else DCTDT,
      DCTADY = if (!("DCTADY" %in% names(pick(everything())))) as.integer(NA) else as.integer(DCTADY),
      LTVISIT = if (!("LTVISIT" %in% names(pick(everything())))) NA_character_ else as.character(LTVISIT)
    )

  reasons_dcts <- c(
    "Adverse Event",
    "Lack of Efficacy",
    "Withdrawal by Subject",
    "Protocol Violation",
    "Death",
    "Other"
  )

  id_vars <- c(
    "STUDYID", "USUBJID", "SITEID", "SUBJID"
  )
  id_vars <- intersect(id_vars, names(adsl))

  adsl <- adsl |>
    mutate(
      TRTEDT_STD  = if ("TRTEDT"  %in% names(pick(everything()))) TRTEDT  else if ("TRT01EDT"  %in% names(pick(everything()))) TRT01EDT  else as.Date(NA),
      TRTEDTM_STD = if ("TRTEDTM" %in% names(pick(everything()))) TRTEDTM else if ("TRT01EDTM" %in% names(pick(everything()))) TRT01EDTM else NA,
      TRTEDY_STD  = if ("TRTEDY"  %in% names(pick(everything()))) as.integer(TRTEDY)  else if ("TRT01EDY"  %in% names(pick(everything()))) as.integer(TRT01EDY)  else as.integer(NA),
      TRTSDT_STD  = if ("TRTSDT"  %in% names(pick(everything()))) TRTSDT  else if ("TRT01SDT"  %in% names(pick(everything()))) TRT01SDT  else as.Date(NA),
      TRTSDTM_STD = if ("TRTSDTM" %in% names(pick(everything()))) TRTSDTM else if ("TRT01SDTM" %in% names(pick(everything()))) TRT01SDTM else NA,
      TRTSDY_STD  = if ("TRTSDY"  %in% names(pick(everything()))) {
        as.integer(TRTSDY)
      } else if ("TRT01SDY"  %in% names(pick(everything()))) {
        as.integer(TRT01SDY)
      } else {
        ifelse(!is.na(TRTSDT_STD), 1L, as.integer(NA))
      }
    )

  eots1 <- adsl |>
    select(any_of(id_vars)) |>
    mutate(
      PARAMCD = "EOTS1STT",
      PARAM = "End of Trt Status for study agent 1",
      AVALC = sample(c("ONGOING", "DISCONTINUED", "COMPLETED"), n(), replace = TRUE),
      DSSCAT = "Study Agent 1"
    )

  eots2 <- adsl |>
    select(any_of(id_vars)) |>
    mutate(
      PARAMCD = "EOTS2STT",
      PARAM = "End of Trt Status for study agent 2",
      AVALC = sample(c("ONGOING", "DISCONTINUED", "COMPLETED"), n(), replace = TRUE),
      DSSCAT = "Study Agent 2"
    )

  dcts1 <- eots1 |>
    filter(AVALC == "DISCONTINUED") |>
    select(any_of(id_vars)) |>
    left_join(
      adsl |>
        transmute(
          across(all_of(id_vars)),
          ASTDT = DCTDT,
          ASTDTC = format(as.Date(DCTDT), "%Y-%m-%d"),
          ASTDY = as.integer(DCTADY)
        ),
      by = id_vars
    ) |>
    mutate(
      PARAMCD = "DCTS1RS",
      PARAM = "Reason for Treatment Discontinuation for study agent 1",
      AVALC = sample(reasons_dcts, n(), replace = TRUE),
      AVALCSP = ifelse(AVALC == "Other", "specify text", NA_character_)
    )

  dcts2 <- eots2 |>
    filter(AVALC == "DISCONTINUED") |>
    select(any_of(id_vars)) |>
    left_join(
      adsl |>
        transmute(
          across(all_of(id_vars)),
          ASTDT = DCTDT,
          ASTDTC = format(as.Date(DCTDT), "%Y-%m-%d"),
          ASTDY = as.integer(DCTADY)
        ),
      by = id_vars
    ) |>
    mutate(
      PARAMCD = "DCTS2RS",
      PARAM = "Reason for Treatment Discontinuation for study agent 2",
      AVALC = sample(reasons_dcts, n(), replace = TRUE),
      AVALCSP = ifelse(AVALC == "Other", "specify text", NA_character_)
    )

  ltv <- adsl |>
    filter(!is.na(LTVISIT) & LTVISIT != "") |>
    select(any_of(id_vars), LTVISIT)

  ltv1 <- ltv |>
    transmute(
      across(all_of(id_vars)),
      PARAMCD = "LTVISTS1",
      PARAM = "Last Treatment Visit for study agent 1",
      AVALC = as.character(LTVISIT)
    )

  ltv2 <- ltv |>
    transmute(
      across(all_of(id_vars)),
      PARAMCD = "LTVISTS2",
      PARAM = "Last Treatment Visit for study agent 2",
      AVALC = as.character(LTVISIT)
    )

  svars <- adsl |>
    transmute(
      across(all_of(id_vars)),
      S1EDT = TRTEDT_STD,
      S1EDTM = TRTEDTM_STD,
      S1EDY = TRTEDY_STD,
      S2EDT = TRTEDT_STD,
      S2EDTM = TRTEDTM_STD,
      S2EDY = TRTEDY_STD,
      S1SDT = TRTSDT_STD,
      S1SDTM = TRTSDTM_STD,
      S1SDY = TRTSDY_STD,
      S2SDT = TRTSDT_STD,
      S2SDTM = TRTSDTM_STD,
      S2SDY = TRTSDY_STD
    )

  add_subj_vars <- c(
    "PPROTFL", "SAFFL", "RANDFL", "SCRNFL", "SCRFFL", "RESCRNFL",
    "FASFL", "ENRLFL", "ENRFL", "AGE", "AGEU", "SEX", "RACE", "COUNTRY",
    "RFICDT", "SITEID", "SUBJID", "TRT01P", "TRT01PN", "TRT01A", "TRT01AN"
  )
  missing_vars <- setdiff(add_subj_vars, names(adsl))
  orig_names <- names(adsl)
  subj <- adsl |>
    mutate(
      !!!setNames(rep(list(NA), length(missing_vars)), missing_vars),
      RFICDT  = if ("RFICDT" %in% names(pick(everything()))) {
        suppressWarnings(format(as.Date(RFICDT), "%Y-%m-%d"))
      } else if ("RFICDTC" %in% names(pick(everything()))) {
        suppressWarnings(format(as.Date(RFICDTC), "%Y-%m-%d"))
      } else if ("DMDTC" %in% names(pick(everything()))) {
        suppressWarnings(format(as.Date(DMDTC), "%Y-%m-%d"))
      } else {
        NA_character_
      },
      PPROTFL = if ("PPROTFL" %in% names(pick(everything()))) {
        if (is.factor(PPROTFL)) as.character(PPROTFL)
        else if (is.logical(PPROTFL)) ifelse(PPROTFL, "Y", "N")
        else if (is.character(PPROTFL)) toupper(PPROTFL)
        else as.character(PPROTFL)
      } else NA_character_,
      RANDFL = if ("RANDFL" %in% names(pick(everything()))) {
        if (is.factor(RANDFL)) as.character(RANDFL)
        else if (is.logical(RANDFL)) ifelse(RANDFL, "Y", "N")
        else if (is.character(RANDFL)) toupper(RANDFL)
        else as.character(RANDFL)
      } else NA_character_,
      ENRLFL = if ("ENRLFL" %in% names(pick(everything()))) {
        if (is.factor(ENRLFL)) as.character(ENRLFL)
        else if (is.character(ENRLFL)) toupper(ENRLFL)
        else as.character(ENRLFL)
      } else NA_character_,
      ENRFL = if ("ENRFL" %in% orig_names) {
        if (is.factor(ENRFL)) as.character(ENRFL)
        else if (is.character(ENRFL)) toupper(ENRFL)
        else as.character(ENRFL)
      } else if ("ENRLFL" %in% orig_names) {
        toupper(as.character(ENRLFL))
      } else NA_character_,
      TRT01PN = if (!"TRT01PN" %in% names(pick(everything()))) as.integer(NA) else as.integer(TRT01PN),
      TRT01AN = if (!"TRT01AN" %in% names(pick(everything()))) as.integer(NA) else as.integer(TRT01AN)
    ) |>
    select(any_of(unique(c(id_vars, add_subj_vars))))

  disp_par <- bind_rows(eots1, dcts1, eots2, dcts2, ltv1, ltv2) |>
    left_join(svars, by = id_vars) |>
    left_join(subj, by = id_vars, suffix = c("", ""))

  additional_labels <- list(
    PARAMCD = "Parameter Code",
    PARAM = "Parameter",
    AVALC    = "Analysis Value (C)",
    AVALCSP = "Specify text for AVALC",
    ASTDT = "Analysis Start Date",
    ASTDTC = "Analysis Start Date",
    ASTDY = "Study Day of Analysis Start Date",
    DSSCAT = "Analysis Subcategory"
  )

  gen <- disp_par |>
    restore_labels(orig_df = adsl, additional_labels = additional_labels, source_dfs = list(adsl))

  return(gen)
}

addisp <- gen_addisp()
