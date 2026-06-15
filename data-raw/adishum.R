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
   # tribble of PARAMCD, PARAM, ISTESTCD
  temp_paramcd <- tibble::tribble(
    ~PARAMCD,  ~PARAM,                                                      ~ISTESTCD,
    "ADABL",    "Binding ADA, Last Result",                                 "ADA_BAB",
    "ADABLT",   "Binding ADA, Titer",                                       "ADA_BAB",
    "ADATRB",   "Treatment-Emergent ADA, Result",                           "ADA_BAB",
    "ADATRBT",  "Treatment-Emergent ADA, Titer",                            "ADA_BAB",
    "ADANTRB",  "Neutralizing ADA, Result",                                 "ADA_NAB",
    "ADANTRBT", "Neutralizing ADA, Titer",                                  "ADA_NAB",
    "ADATRI",   "Treatment-Induced ADA, Result",                            "ADA_BAB",
    "ADATRIPT", "Treatment-Induced ADA, Titer",                             "ADA_BAB",
    "ADATRE",   "Treatment-Enhanced ADA, Result",                           "ADA_BAB",
    "ADATREPT", "Treatment-Enhanced ADA, Titer",                            "ADA_BAB",
    "ADANTRE",  "Treatment-Enhanced Neutralizing ADA, Result",              "ADA_NAB",
    "ADAPSP",   "ADA Persistent Positive, Result",                          "ADA_BAB",
    "ADATSP",   "Neutralizing ADA Persistent Positive",                     "ADA_NAB",
    "ADAUND",   "ADA Undetermined, Result",                                 "ADA_BAB",
    "TMOSADAW", "Time to Onset of ADA (weeks)",                             "ADA_BAB",
    "ADADURW",  "Duration of ADA Response (weeks)",                         "ADA_BAB",
    "PSPDURW",  "Duration of Persistent Positive ADA (weeks)",              "ADA_BAB",
    "NABPOS",   "Neutralizing ADA Positive",                                "ADA_NAB",
    "NABNEG",   "Neutralizing ADA Negative",                                "ADA_NAB"
  ) |> 
    arrange(PARAMCD)

  # initialize columns
  main_tbl <- main_tbl |> 
    mutate(
      PARAMCD = NA_character_,
      PARAM   = NA_character_
    )

  for (val in c("ADA_BAB", "ADA_NAB")) {
    idx  <- which(as.character(main_tbl$ISTESTCD) == val)     # rows to fill
    pool <- temp_paramcd |> filter(ISTESTCD == val)          # candidate PARAM rows

    sel <- sample(
      seq_along(pool$PARAMCD),
      size = length(idx),
      replace = TRUE
    )

    main_tbl$PARAMCD[idx] <- pool$PARAMCD[sel]
    main_tbl$PARAM[idx]   <- pool$PARAM[sel]
  }

  return(main_tbl)

}

# add AVAL, AVALC and AVALCAT1 columns
add_adishum_col_funcs$aval_avalc_avalcat1 <- function(main_tbl) {
  # Group B parameters (Subject Summary) handled later
  group_b_paramcd <- c("TMOSADAW", "ADADURW", "PSPDURW", "NABPOS", "NABNEG")

  # Row-wise small noise for numeric AVAL
  noise <- stats::runif(n = nrow(main_tbl), min = -0.05, max = 0.05)

  main_tbl <- main_tbl |>
    dplyr::mutate(
      # Numeric AVAL for ISSTRESN and not Group B, rounded to 2 decimals (spec: 2 significant digits)
      AVAL = dplyr::case_when(
        !is.na(ISSTRESN) & !PARAMCD %in% group_b_paramcd ~ round(ISSTRESN + noise, 2),
        TRUE ~ NA_real_
      ),
      # AVALC mirrors numeric AVAL as character; else keep character ISSTRESC (upper-case)
      AVALC = dplyr::case_when(
        !is.na(AVAL) & !PARAMCD %in% group_b_paramcd ~ as.character(AVAL),
        is.na(ISSTRESN) & !is.na(ISSTRESC) & !PARAMCD %in% group_b_paramcd ~ toupper(ISSTRESC),
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::mutate(
      # AVALCAT1 per spec (Group B set later)
      AVALCAT1 = dplyr::case_when(
        !PARAMCD %in% group_b_paramcd &
          ( (!is.na(AVAL) & AVAL <= 0) |
            (!is.na(AVALC) & AVALC == "NEGATIVE") |
            (!is.na(AVALC) & grepl("^\\s*<", AVALC)) ) ~ "Negative / BLQ",
        !PARAMCD %in% group_b_paramcd & !is.na(AVAL) & AVAL > 0 & AVAL <= 2 ~ "Low Positive",
        !PARAMCD %in% group_b_paramcd & !is.na(AVAL) & AVAL > 2 ~ "High Positive",
        TRUE ~ NA_character_
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
  # add avisitn
  main_tbl <- main_tbl |>
    dplyr::arrange(USUBJID, ISDTC) |>
    dplyr::group_by(USUBJID) |>
    dplyr::mutate(AVISITN = dplyr::row_number()) |>
    dplyr::ungroup()

  # add adt, ady, ablfl
  main_tbl <- main_tbl |>
    dplyr::mutate(
      ADT = as.Date(substr(ISDTC, 1, 10)),
      ADY = dplyr::if_else(
        ADT >= TRTSDT,
        as.integer(ADT - TRTSDT + 1L),
        as.integer(ADT - TRTSDT)
      ),
      ABLFL = dplyr::if_else(
        AVISITN == 1, "Y", NA_character_
      )
    )

  main_tbl <- main_tbl |>
    dplyr::mutate(
      AVISIT = dplyr::case_when(
        ABLFL == "Y"  ~ "Baseline",
        AVISITN == 1  ~ "Week 2",
        AVISITN == 2  ~ "Week 4",
        AVISITN == 3  ~ "Week 8",
        AVISITN == 4  ~ "Week 12",
        AVISITN == 5  ~ "Week 16",
        AVISITN == 6  ~ "Week 20",
        AVISITN >= 7  ~ "Week 24",
        TRUE          ~ NA_character_
      )
    )

  return(main_tbl)
}

# add ANL01FL and ANL02FL columns
add_adishum_col_funcs$anl01fl_anl02fl <- function(main_tbl) {

  # target study days per visit number (proxy: every visit ~2 weeks apart)
  # covers up to AVISITN = 20 to be safe since AVISITN is row_number() per subject
  visitdy_map <- setNames(
    as.numeric(seq(15, by = 14, length.out = 19)),  # day 15, 29, 43 ... up to visit 20
    as.character(2:20)                               # visit 2 to 20 (visit 1 = baseline)
  )

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
      .tgt = dplyr::recode(as.character(AVISITN), !!!visitdy_map, .default = NA_real_)
    ) |>
    dplyr::group_by(USUBJID, PARAMCD, AVISITN) |>
    dplyr::mutate(
      .eligible = ANL01FL == "Y" & AVISITN > 1 & !is.na(.tgt),
      .min_diff = if (any(.eligible, na.rm = TRUE))
                    min(abs(ADY[.eligible] - .tgt), na.rm = TRUE)
                  else
                    NA_real_,
      ANL02FL = dplyr::if_else(
        .eligible & abs(ADY - .tgt) == .min_diff,
        "Y", NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.tgt, -.eligible, -.min_diff)

  return(main_tbl)
}


# core function ------------------------------
# Generate adishum dataset
gen_adishum <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # get source data
  sdtm_tbl <- pharmaversesdtm::is_ada |> 
    left_join(
      pharmaverseadamjnj::adsl,
      by = "USUBJID"
    ) |> 
    select(
      USUBJID,
      SITEID, 
      SUBJID,
      AGE,
      SEX,
      RACE,
      IMFL,
      SAFFL,
      ISBDAGNT,
      TRTSDT, 
      TRTSDTM,
      TRT01P, 
      TRT01PN, 
      TRT01A, 
      TRT01AN,
      ISSTRESN,
      ISSTRESC,
      ISTESTCD,
      ISDTC
    )
  
  sdtm_tbl <- sdtm_tbl |> 
    left_join(
      pharmaverseadam::adab |>
        filter(ABLFL == "Y") |>
        select(USUBJID, ADTM) |>
        distinct(USUBJID, .keep_all = TRUE),
      by = "USUBJID"
    )

  # add PARQUAL column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_adishum_col_funcs$parqual()

  # add PARAM and PARAMCD column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_adishum_col_funcs$param_n_paramcd()

  # add AVAL, AVALC and AVALCAT1 column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_adishum_col_funcs$aval_avalc_avalcat1()

  # add IMEVFL column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_adishum_col_funcs$imevfl()

  # add ADT, ADY, ABLFL, AVISIT and AVISITN column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_adishum_col_funcs$adt_ady_ablfl_avisit_avisitn()

  # add ANL01FL and ANL02FL column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_adishum_col_funcs$anl01fl_anl02fl()

}


# run function ------------------------------
# generate dataset
adishum <- gen_adishum()
