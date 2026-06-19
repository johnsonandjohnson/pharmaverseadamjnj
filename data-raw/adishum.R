# generate adishum dataset
# setup section ------------------------------

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(formatters)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# helper functions ----------------------------
# functions to keep the code readable and tidy

add_adishum_col_funcs <- list()

# function for adding sample visits
add_adishum_col_funcs$sample_visits <- function() {
  visit_map <- tibble::tribble(
    ~VISITNUM, ~AVISITN, ~AVISIT, ~VISITDY,
    1L, 1L, "Screening", -14,
    2L, 2L, "Day 1", 1,
    3L, 3L, "Baseline", 1,
    4L, 4L, "Week 2", 15,
    5L, 5L, "Week 4", 29,
    6L, 6L, "Week 6", 43,
    7L, 7L, "Week 8", 57,
    8L, 8L, "Week 10", 71,
    9L, 9L, "Week 12", 85,
    10L, 10L, "Week 16", 113,
    11L, 11L, "Week 20", 141,
    12L, 12L, "Week 24", 169,
    13L, 13L, "Week 32", 225,
    14L, 14L, "Week 40", 281,
    15L, 15L, "Week 48", 337
  )

  return(visit_map)
}

# function for adding sample paramcd
add_adishum_col_funcs$sample_paramcd <- function() {
  # tribble of PARAMCD, PARAM, ISTESTCD
  sample_paramcd <- tibble::tribble(
    ~PARAMCD, ~PARAM, ~ISTESTCD,
    "ADABL", "Binding ADA, Last Result", "ADA_BAB",
    "ADABLT", "Binding ADA, Titer", "ADA_BAB",
    "ADATRB", "Treatment-Emergent ADA, Result", "ADA_BAB",
    "ADATRBT", "Treatment-Emergent ADA, Titer", "ADA_BAB",
    "ADANTRB", "Neutralizing ADA, Result", "ADA_NAB",
    "ADANTRBT", "Neutralizing ADA, Titer", "ADA_NAB",
    "ADATRI", "Treatment-Induced ADA, Result", "ADA_BAB",
    "ADATRIPT", "Treatment-Induced ADA, Titer", "ADA_BAB",
    "ADATRE", "Treatment-Enhanced ADA, Result", "ADA_BAB",
    "ADATREPT", "Treatment-Enhanced ADA, Titer", "ADA_BAB",
    "ADANTRE", "Treatment-Enhanced Neutralizing ADA, Result", "ADA_NAB",
    "ADAPSP", "ADA Persistent Positive, Result", "ADA_BAB",
    "ADATSP", "Neutralizing ADA Persistent Positive", "ADA_NAB",
    "ADAUND", "ADA Undetermined, Result", "ADA_BAB",
    "TMOSADAW", "Time to Onset of ADA (weeks)", NA_character_,
    "ADADURW", "Duration of ADA Response (weeks)", NA_character_,
    "PSPDURW", "Duration of Persistent Positive ADA (weeks)", NA_character_,
    "NABPOS", "Neutralizing ADA Positive", NA_character_,
    "NABNEG", "Neutralizing ADA Negative", NA_character_
  ) |>
    arrange(PARAMCD)

  return(sample_paramcd)
}

# function for adding PARQUAL col
add_adishum_col_funcs$parqual <- function(main_tbl) {
  main_tbl <- main_tbl |>
    mutate(
      PARQUAL = if_else(
        is.na(ISBDAGNT),
        NA_character_,
        ISBDAGNT |>
          stringi::stri_trim_both() |>
          stringi::stri_trans_toupper()
      )
    )

  return(main_tbl)
}

# function for adding PARAM and PARAMCD col
add_adishum_col_funcs$param_n_paramcd <- function(main_tbl) {
  sample_paramcd <- add_adishum_col_funcs$sample_paramcd()
  # initialize columns
  main_tbl <- main_tbl |>
    mutate(
      PARAMCD = NA_character_,
      PARAM = NA_character_
    )

  # for ADA_BAB and ADA_NAB in ISTESTCD
  for (val in c("ADA_BAB", "ADA_NAB")) {
    idx <- which(as.character(main_tbl$ISTESTCD) == val) # rows to fill
    pool <- sample_paramcd |> filter(ISTESTCD == val) # candidate PARAM rows

    sel <- sample(
      seq_along(pool$PARAMCD),
      size = length(idx),
      replace = TRUE
    )

    main_tbl$PARAMCD[idx] <- pool$PARAMCD[sel]
    main_tbl$PARAM[idx] <- pool$PARAM[sel]
  }

  # For NA values in ISTESTCD
  na_paramcd <- sample_paramcd |>
    dplyr::filter(is.na(ISTESTCD))

  if (nrow(na_paramcd) > 0) {
    sel_na <- sample(
      x = seq_along(main_tbl$USUBJID),
      size = floor(nrow(main_tbl) * 0.2),
      replace = FALSE
    )

    na_paramcd <- na_paramcd[
      rep(
        na_paramcd$PARAMCD |> seq_along(),
        (nrow(main_tbl) / nrow(na_paramcd)) |>
          ceiling()
      ),
    ]

    na_rows_tbl <- main_tbl[sel_na, ]

    na_rows_tbl$ISTESTCD <- NA_character_
    na_rows_tbl$PARAMCD <- na_paramcd$PARAMCD[sel_na]
    na_rows_tbl$PARAM <- na_paramcd$PARAM[sel_na]

    main_tbl <- main_tbl |>
      dplyr::bind_rows(na_rows_tbl)
  }

  return(main_tbl)
}

# add AVAL, AVALC and AVALCAT1 columns
add_adishum_col_funcs$aval_avalc_avalcat1 <- function(main_tbl) {
  # Group B parameters (Subject Summary) handled later
  group_b_paramcd <- c("TMOSADAW", "ADADURW", "PSPDURW", "NABPOS", "NABNEG")

  # initialize columns
  main_tbl <- main_tbl |>
    dplyr::mutate(
      AVAL = NA_real_,
      AVALC = NA_character_,
      AVALCAT1 = NA_character_
    )

  # override AVAL for ADATREPT to ensure all 4 titer AVALCAT1 tiers are covered
  adatrept_idx <- which(main_tbl$PARAMCD == "ADATREPT")

  if (length(adatrept_idx) > 0) {
    main_tbl$AVAL[adatrept_idx] <- sample(
      c(1:9, 91:99, 991:999, 1001:1009),
      length(adatrept_idx),
      replace = TRUE
    )
  }

  # Row-wise small noise for numeric AVAL
  noise <- stats::runif(n = nrow(main_tbl), min = -0.05, max = 0.05)

  main_tbl <- main_tbl |>
    dplyr::mutate(
      AVAL = dplyr::case_when(
        PARAMCD == "ADATREPT" ~ AVAL,
        !is.na(ISSTRESN) & !PARAMCD %in% group_b_paramcd ~ round(ISSTRESN + noise, 2)
      ),
      AVALC = dplyr::case_when(
        !is.na(AVAL) & !PARAMCD %in% group_b_paramcd ~ as.character(AVAL),
        is.na(ISSTRESN) & !is.na(ISSTRESC) & !PARAMCD %in% group_b_paramcd ~ toupper(ISSTRESC)
      )
    ) |>
    dplyr::mutate(
      AVALCAT1 = dplyr::case_when(
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL < 10 ~ "<10",
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL < 100 ~ "10 to < 100",
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL < 1000 ~ "100 to <1000",
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL >= 1000 ~ ">= 1000",
        !PARAMCD %in% group_b_paramcd &
          ((!is.na(AVAL) & AVAL <= 0) |
            (!is.na(AVALC) & AVALC == "NEGATIVE") |
            (!is.na(AVALC) & grepl("^\\s*<", AVALC))) ~ "Negative / BLQ",
        !PARAMCD %in% group_b_paramcd & !is.na(AVAL) & AVAL > 0 & AVAL <= 2 ~ "Low Positive",
        !PARAMCD %in% group_b_paramcd & !is.na(AVAL) & AVAL > 2 ~ "High Positive"
      )
    )

  return(main_tbl)
}


# add IMEVFL column
add_adishum_col_funcs$imevfl <- function(main_tbl) {
  # Subject-level eligibility from IMFL
  subj <- main_tbl |>
    dplyr::group_by(USUBJID) |>
    dplyr::summarise(
      eligible = any(as.character(IMFL) == "Y", na.rm = TRUE),
      .groups = "drop"
    )

  # Sample ~70% of eligible subjects as Y
  eligible_ids <- subj |>
    dplyr::filter(eligible) |>
    dplyr::pull(USUBJID)

  n_y <- floor(length(eligible_ids) * 0.70)
  sampled_ids <- if (n_y > 0) {
    sample(eligible_ids, size = n_y)
  } else {
    character(0)
  }

  # Build subject-level flag
  subj_flag <- subj |>
    dplyr::mutate(
      IMEVFL = ifelse(USUBJID %in% sampled_ids, "Y", "N")
    ) |>
    dplyr::select(USUBJID, IMEVFL)

  # Attach to main table, adding only IMEVFL
  main_tbl <- main_tbl |>
    dplyr::left_join(
      subj_flag,
      by = "USUBJID"
    )

  return(main_tbl)
}

# add AVISIT and AVISITN columns
add_adishum_col_funcs$adt_ady_ablfl_avisit_avisitn <- function(main_tbl) {
  main_tbl <- main_tbl |>
    dplyr::mutate(
      AVISITN = VISITNUM,
      AVISIT = stringr::str_to_title(VISIT),
      ADT = as.Date(substr(ISDTC, 1, 10))
    )

  # added ADY
  main_tbl <- main_tbl |>
    dplyr::mutate(
      ADY = dplyr::if_else(
        ADT >= TRTSDT,
        as.integer(ADT - TRTSDT + 1L),
        as.integer(ADT - TRTSDT)
      )
    )

  # added ABLFL
  main_tbl <- main_tbl |>
    dplyr::group_by(USUBJID, PARAMCD) |>
    dplyr::mutate(
      ABLFL = dplyr::if_else(
        VISITNUM == min(VISITNUM, na.rm = TRUE),
        "Y",
        NA_character_
      )
    ) |>
    dplyr::ungroup()

  return(main_tbl)
}

# add ANL01FL and ANL02FL columns
add_adishum_col_funcs$anl01fl_anl02fl <- function(main_tbl) {
  # get sample visit mapping
  vm <- add_adishum_col_funcs$sample_visits()
  visitdy_map <- setNames(as.numeric(vm$VISITDY), as.character(vm$AVISITN))

  main_tbl <- main_tbl |>
    # ensure consistent row order before any row_number() logic
    dplyr::arrange(USUBJID, PARAMCD, ADT, ADTM) |>
    # ANL01FL: last non-missing record per subject + parameter + date
    dplyr::group_by(USUBJID, PARAMCD, ADT) |>
    dplyr::mutate(
      ANL01FL = dplyr::if_else(
        (!is.na(AVAL) | !is.na(AVALC)) & dplyr::row_number() == dplyr::n(),
        "Y",
        NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    # ANL02FL: one record per scheduled post-baseline visit,
    #          closest actual day (ADY) to the visit target day
    dplyr::mutate(
      .tgt = dplyr::recode(
        as.character(AVISITN),
        !!!visitdy_map,
        .default = NA_real_
      )
    ) |>
    dplyr::group_by(USUBJID, PARAMCD, AVISITN) |>
    dplyr::mutate(
      .eligible = ANL01FL == "Y" & AVISITN > 1 & !is.na(.tgt),
      .min_diff = if (any(.eligible, na.rm = TRUE)) {
        min(abs(ADY[.eligible] - .tgt), na.rm = TRUE)
      } else {
        NA_real_
      },
      ANL02FL = dplyr::if_else(
        .eligible & abs(ADY - .tgt) == .min_diff,
        "Y",
        NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.tgt, -.eligible, -.min_diff)

  return(main_tbl)
}

# add ANL03FL, ANL04FL, ANL05FL, ANL06FL, ANL07FL, ANL08FL, ANL09FL and ANL10FL columns
add_adishum_col_funcs$anl03fl_to_anl10fl <- function(main_tbl) {
  group_b_paramcd <- c("TMOSADAW", "ADADURW", "PSPDURW", "NABPOS", "NABNEG")

  # subject-level flags — built once, joined back
  subj_flags <- main_tbl |>
    dplyr::filter(!is.na(TRTSDT)) |>
    dplyr::distinct(USUBJID)

  subj_flags <- subj_flags |>
    dplyr::mutate(
      # random 12% sampling for the column
      ANL03FL = dplyr::if_else(
        dplyr::row_number() %in% sample(dplyr::n(), floor(dplyr::n() * 0.12)),
        "Y",
        NA_character_
      ),
      # random 10% sampling for the column
      ANL07FL = dplyr::if_else(
        dplyr::row_number() %in% sample(dplyr::n(), floor(dplyr::n() * 0.10)),
        "Y",
        NA_character_
      )
    )

  # all derived flags
  anl03_y_idx <- which(subj_flags$ANL03FL == "Y")
  anl07_y_idx <- which(subj_flags$ANL07FL == "Y")

  # take percetages of rows in random order.
  subj_flags <- subj_flags |>
    dplyr::mutate(
      ANL04FL = dplyr::if_else(
        row_number() %in% sample(anl03_y_idx, floor(length(anl03_y_idx) * 0.35)),
        "Y",
        NA_character_
      ),
      ANL05FL = dplyr::if_else(
        row_number() %in% sample(anl03_y_idx, floor(length(anl03_y_idx) * 0.25)),
        "Y",
        NA_character_
      ),
      ANL06FL = dplyr::if_else(
        row_number() %in% sample(anl03_y_idx, floor(length(anl03_y_idx) * 0.15)),
        "Y",
        NA_character_
      ),
      ANL08FL = dplyr::if_else(
        row_number() %in% sample(anl07_y_idx, floor(length(anl07_y_idx) * 0.35)),
        "Y",
        NA_character_
      ),
      ANL09FL = dplyr::if_else(
        row_number() %in% sample(anl07_y_idx, floor(length(anl07_y_idx) * 0.25)),
        "Y",
        NA_character_
      ),
      ANL10FL = dplyr::if_else(
        row_number() %in% sample(anl07_y_idx, floor(length(anl07_y_idx) * 0.15)),
        "Y",
        NA_character_
      )
    )

  # join back to main table
  main_tbl <- main_tbl |>
    dplyr::left_join(subj_flags, by = "USUBJID") |>
    # group B (subject summary) rows get NA for all reaction flags
    dplyr::mutate(
      dplyr::across(
        ANL03FL:ANL10FL,
        ~ dplyr::if_else(
          PARAMCD %in% group_b_paramcd,
          NA_character_,
          .x
        )
      )
    )

  return(main_tbl)
}

# ensure required AVALC values are present per PARAMCD spec
add_adishum_col_funcs$pp_ensure_required_avalc <- function(
  main_tbl,
  min_rows = 5
) {
  y_paramcds <- c(
    "ADABL",
    "ADATRB",
    "ADANTRB",
    "ADATRI",
    "ADATRE",
    "ADANTRE",
    "ADAPSP",
    "ADATSP",
    "ADAUND",
    "NABPOS",
    "NABNEG"
  )
  adatrept_cats <- c("<10", "10 to < 100", "100 to <1000", ">= 1000")

  # For each PARAMCD requiring AVALC="Y", force min_rows random rows to "Y"
  main_tbl <- main_tbl |>
    dplyr::group_by(USUBJID, PARAMCD) |>
    dplyr::mutate(
      AVALC = dplyr::if_else(
        {
          PARAMCD %in%
            y_paramcds &
            dplyr::row_number() %in%
              sample(dplyr::n(), min(min_rows, dplyr::n()))
        },
        "Y",
        AVALC
      )
    ) |>
    dplyr::ungroup()

  return(main_tbl)
}

# core function ------------------------------
# Generate adishum dataset
gen_adishum <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # get available data --------------
  # get source data
  raw <- pharmaversesdtm::is_ada

  # Variables to keep from IS source
  is_ada_vars <- c(
    "USUBJID",
    "ISSEQ",
    "ISTESTCD",
    "ISTEST",
    "ISBDAGNT",
    "ISCAT",
    "ISSTRESC",
    "ISSTRESN",
    "ISNAM",
    "ISSPEC",
    "VISITNUM",
    "VISIT",
    "ISDTC"
  )

  # Variables to keep exclusively from ADSL
  adsl_vars <- c(
    "STUDYID",
    "SUBJID",
    "SITEID",
    "AGE",
    "SEX",
    "RACE",
    "IMFL",
    "SAFFL",
    "TRTSDT",
    "TRTSDTM",
    "TRT01P",
    "TRT01PN",
    "TRT01A",
    "TRT01AN"
  )

  adsl_subset <- pharmaverseadamjnj::adsl |>
    dplyr::select(
      USUBJID,
      dplyr::all_of(adsl_vars)
    )

  gen <- raw |>
    dplyr::select(
      dplyr::all_of(is_ada_vars)
    ) |>
    dplyr::left_join(
      adsl_subset,
      by = "USUBJID"
    )

  gen <- gen |>
    dplyr::mutate(
      ADTM = admiral::convert_dtc_to_dtm(ISDTC)
    )

  # add new variables --------------
  # add PARQUAL columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$parqual()

  # add PARAM and PARAMCD columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$param_n_paramcd()

  # add AVAL, AVALC and AVALCAT1 columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$aval_avalc_avalcat1()

  # add IMEVFL columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$imevfl()

  # add ADT, ADY, ABLFL, AVISIT and AVISITN columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$adt_ady_ablfl_avisit_avisitn()

  # add ANL01FL and ANL02FL columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$anl01fl_anl02fl()

  # add ANL03FL, ANL04FL, ANL05FL, ANL06FL, ANL07FL, ANL08FL, ANL09FL and ANL10FL columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$anl03fl_to_anl10fl()

  # ensure required AVALC values are present per PARAMCD spec
  gen <- gen |>
    add_adishum_col_funcs$pp_ensure_required_avalc()

  # Handle NA values and convert characters to factors
  gen <- gen |>
    df_na(
      char_as_factor = TRUE
    )

  # Define additional labels for new variables not in source dataset
  # Mrinal's labels take precedence — these override all other sources
  additional_labels <- list(
    IMFL = "Immunogenicity Population Flag",
    PARQUAL = "Parameter Qualifier",
    PARAMCD = "Parameter Code",
    AVAL = "Analysis Value",
    AVALC = "Analysis Value (C)",
    AVALCAT1 = "Analysis Value Category 1",
    IMEVFL = "IS Evaluable Population Flag",
    ANL01FL = "Analysis Flag 01",
    ANL02FL = "Analysis Flag 02",
    ANL03FL = "Analysis Flag 03-Inf SR",
    ANL04FL = "Analysis Flag 04-Severe Inf SR",
    ANL05FL = "Analysis Flag 05-Serious Inf SR",
    ANL06FL = "Analysis Flag 06-Infus SR - DC",
    ANL07FL = "Analysis Flag 07-Inj SR",
    ANL08FL = "Analysis Flag 08-Severe Inj SR",
    ANL09FL = "Analysis Flag 09-Serious Inj SR",
    ANL10FL = "Analysis Flag 10-Inject SR - DC",
    AVISIT = "Analysis Visit",
    AVISITN = "Analysis Visit (N)",
    PARAM = "Parameter Description",
    ADT = "Analysis Date",
    ADY = "Analysis Relative Day",
    ABLFL = "Baseline Record Flag",
    ADTM = "Analysis Datetime",
    SUBJID = "Subject Identifier for the Study",
    SITEID = "Study Site Identifier",
    AGE = "Age",
    SEX = "Sex",
    RACE = "Race",
    SAFFL = "Safety Population Flag",
    TRTSDT = "Date of First Exposure to Treatment",
    TRTSDTM = "Datetime of First Exposure to Treatment",
    TRT01P = "Planned Treatment for Period 01",
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01A = "Actual Treatment for Period 01",
    TRT01AN = "Actual Treatment for Period 01 (N)"
  )

  # Restore labels
  gen <- restore_labels(
    df = gen,
    orig_df = raw,
    additional_labels = additional_labels
  )

  # Return Gen
  return(gen)
}

# run function ------------------------------
# generate dataset
adishum <- gen_adishum(seed = 123)
