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
  tibble::tribble(
    ~VISITNUM , ~AVISITN , ~AVISIT , ~VISITDY ,
     1L       , 1L       , "Day 1" ,      -14 ,
     2L       , 2L       , "Day 2" ,        1 ,
     3L       , 3L       , "Day 3" ,        1 ,
     4L       , 1L       , "Day 1" ,       15 ,
     5L       , 2L       , "Day 2" ,       29 ,
     6L       , 3L       , "Day 3" ,       43 ,
     7L       , 1L       , "Day 1" ,       57 ,
     8L       , 2L       , "Day 2" ,       71 ,
    12L       , 3L       , "Day 3" ,      169
  )
}

# function for adding sample paramcd
add_adishum_col_funcs$sample_paramcd <- function() {
  # tribble of PARAMCD, PARAM, ISTESTCD
  sample_param_tbl <- tibble::tribble(
    ~VISITNUM , ~ISTESTCD     , ~PARCAT1          , ~PARAMCD   , ~PARAM                                                 ,
            1 , "ADA_BAB"     , "Collection"      , "SCRRSLT"  , "Screening Result"                                     ,
            2 , "ADA_BAB"     , "Collection"      , "CNRRSLT"  , "Confirmatory Result"                                  ,
            3 , "ADA_BAB"     , "Collection"      , "TITER"    , "Titer"                                                ,
            4 , "ADA_BAB"     , "Collection"      , "SUADAPT"  , "Subject Peak Titer"                                   ,
            5 , "ADA_BAB"     , "Collection"      , "SUADAST"  , "Subject ADA Status"                                   ,
            6 , "ADA_NAB"     , "Collection"      , "NABSCR"   , "Neutralizing Screening Result"                        ,
            7 , "ADA_NAB"     , "Collection"      , "NABCNR"   , "Neutralizing Confirmatory Result"                     ,
            8 , "ADA_NAB"     , "Collection"      , "SUNABST"  , "Subject NAB Status"                                   ,
            0 , NA_character_ , "Subject Summary" , "ADABL"    , "Baseline ADA Positive"                                ,
            0 , NA_character_ , "Subject Summary" , "ADABLT"   , "Baseline ADA Positive Titers"                         ,
            0 , NA_character_ , "Subject Summary" , "ADATRB"   , "Treatment-boosted ADA Positive"                       ,
            0 , NA_character_ , "Subject Summary" , "ADATRBT"  , "Treatment-boosted ADA Positive Titers"                ,
            0 , NA_character_ , "Subject Summary" , "ADANTRB"  , "Not Treatment-boosted ADA Positive"                   ,
            0 , NA_character_ , "Subject Summary" , "ADANTRBT" , "Not Treatment-boosted ADA Positive Titers"            ,
            0 , NA_character_ , "Subject Summary" , "ADATRI"   , "Treatment-induced ADA Positive"                       ,
            0 , NA_character_ , "Subject Summary" , "ADATRIPT" , "Treatment-induced ADA Positive Peak Titers"           ,
            0 , NA_character_ , "Subject Summary" , "ADATRE"   , "Treatment-emergent ADA Positive"                      ,
            0 , NA_character_ , "Subject Summary" , "ADATREPT" , "Treatment-emergent ADA Positive Peak Titers"          ,
            0 , NA_character_ , "Subject Summary" , "ADANTRE"  , "Treatment-emergent ADA Negative"                      ,
            0 , NA_character_ , "Subject Summary" , "ADAPSP"   , "Persistent ADA Response"                              ,
            0 , NA_character_ , "Subject Summary" , "ADATSP"   , "Transient ADA Response"                               ,
            0 , NA_character_ , "Subject Summary" , "ADAUND"   , "Undetermined ADA Response"                            ,
            0 , NA_character_ , "Subject Summary" , "TMOSADAW" , "Time to onset of treatment-induced ADA (weeks)"       ,
            0 , NA_character_ , "Subject Summary" , "ADADURW"  , "Duration of treatment-induced ADA (weeks)"            ,
            0 , NA_character_ , "Subject Summary" , "PSPDURW"  , "Duration of persistent treatment-induced ADA (weeks)" ,
            0 , NA_character_ , "Subject Summary" , "NABPOS"   , "NAB Positive"                                         ,
            0 , NA_character_ , "Subject Summary" , "NABNEG"   , "NAB Negative"                                         ,
  ) |>
    arrange(
      VISITNUM,
      PARAMCD
    )

  return(sample_param_tbl)
}

# sync adsl vars by dropping and re-joining from adsl
add_adishum_col_funcs$sync_adsl_vars <- function(
  main_tbl,
  adsl_vars,
  adsl_dataset
) {
  adsl_subset <- adsl_dataset |>
    dplyr::select(USUBJID, dplyr::all_of(adsl_vars))

  main_tbl <- main_tbl |>
    dplyr::select(-dplyr::any_of(c(adsl_vars, "TRTA"))) |>
    dplyr::left_join(adsl_subset, by = "USUBJID") |>
    dplyr::mutate(TRTA = TRT01A)

  return(main_tbl)
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

  # initialize PARAMCD / PARAM
  main_tbl <- main_tbl |>
    dplyr::mutate(
      PARAMCD = NA_character_,
      PARAM = NA_character_
    )

  # all columns from sample_paramcd that will be overwritten (excl. VISITNUM)
  paramcd_cols <- c(
    "VISITNUM",
    "ISTESTCD",
    "PARCAT1",
    "PARAMCD",
    "PARAM"
  )

  # for Collection: ----------------------------------
  # Collection lookup: SEQ 1-8, one row per SEQ
  collection_seq <- sample_paramcd |>
    dplyr::filter(
      stringi::stri_trans_tolower(PARCAT1) == "collection"
    ) |>
    dplyr::select(VISITNUM, dplyr::all_of(paramcd_cols))

  # --- Step 1: artificial rows for 12 sampled subjects x visitnum 1:8 ---
  # Done BEFORE PARAMCD assignment so the join below covers them too.
  # Clone one existing row per subject per visit, overwrite VISITNUM + VISIT.
  select_subj_n <- 180L
  art_subjects <- main_tbl$USUBJID |>
    unique() |>
    sample(select_subj_n)

  art_rows <- lapply(art_subjects, function(subj) {
    donor <- main_tbl[main_tbl$USUBJID == subj, ]
    out <- donor[sample(nrow(donor), 8L, replace = TRUE), ]
    rownames(out) <- NULL
    out$VISITNUM <- 1:8
    return(out)
  }) |>
    dplyr::bind_rows()

  # bind and keep 1 row per subject x visitnum; existing rows take priority
  main_tbl <- dplyr::bind_rows(main_tbl, art_rows) |>
    dplyr::distinct(
      USUBJID,
      VISITNUM,
      .keep_all = TRUE
    ) |>
    filter(
      VISITNUM <= 9 |
        VISITNUM == 12
    )

  # --- Step 2: Collection rows ---
  # Join on VISITNUM = VISITNUM and overwrite ALL paramcd_cols from collection_seq.
  # Rows with VISITNUM outside 1:8 find no match and keep NA.
  main_tbl <- main_tbl |>
    dplyr::left_join(
      collection_seq,
      by = c("VISITNUM"),
      suffix = c("_old", "")
    ) |>
    dplyr::select(-dplyr::ends_with("_old"))

  # for Subject Summary: ----------------------------------
  # Subject Summary rows: one row per subject × PARAMCD combination
  subj_summary_pool <- sample_paramcd |>
    dplyr::filter(
      stringi::stri_trans_tolower(PARCAT1) == "subject summary"
    )

  # columns whose values come from sample_paramcd (excl. SEQ)
  paramcd_cols <- c(
    "VISITNUM",
    "ISTESTCD",
    "PARCAT1",
    "PARAMCD",
    "PARAM"
  )

  # columns to sample from main_tbl (excl. subject, visit, and paramcd_cols)
  visit_cols <- c("VISIT", "AVISITN", "AVISIT")
  treatment_cols <- c("TRTSDT", "TRTSDTM", "TRT01P", "TRT01PN", "TRT01A", "TRT01AN", "TRTA")
  fixed_cols <- c("USUBJID", visit_cols, paramcd_cols, treatment_cols)
  sample_cols <- setdiff(names(main_tbl), fixed_cols)

  n_ss <- nrow(subj_summary_pool) # total Subject Summary PARAMCDs

  # build per-subject blocks vectorized via lapply, then bind once
  subjects <- unique(main_tbl$USUBJID)

  ss_rows <- lapply(subjects, function(subj) {
    n_pick <- sample(2L:n_ss, 1L) # random count [2, n_ss]
    picked <- subj_summary_pool[sample(n_ss, n_pick), ] # unique PARAMCDs

    # one donor row from main_tbl for this subject to copy non-key cols from
    donor <- main_tbl |>
      filter(USUBJID == subj)

    # build skeleton: n_pick rows with subject fixed
    out <- donor[sample(nrow(donor), n_pick, replace = TRUE), ]
    rownames(out) <- NULL

    # overwrite paramcd_cols from the picked Subject Summary rows
    for (col in paramcd_cols) {
      out[[col]] <- picked[[col]]
    }

    # blank out visit columns (Subject Summary has no visit context)
    for (col in intersect(visit_cols, names(out))) {
      out[[col]] <- NA
    }

    # sample each remaining column from all non-NA values in main_tbl
    for (col in sample_cols) {
      pool_vals <- main_tbl[[col]][!is.na(main_tbl[[col]])]
      if (length(pool_vals) > 0) {
        out[[col]] <- sample(pool_vals, n_pick, replace = TRUE)
      }
    }

    out
  })

  main_tbl <- main_tbl |>
    dplyr::bind_rows(ss_rows)

  return(main_tbl)
}

# add AVAL, AVALC and AVALCAT1 columns
add_adishum_col_funcs$aval_avalc_avalcat1 <- function(main_tbl) {
  group_b_paramcd <- c("TMOSADAW", "ADADURW", "PSPDURW", "NABPOS", "NABNEG")
  titer_paramcds <- c("ADABLT", "ADATRBT", "ADANTRBT", "ADATRIPT", "ADATREPT")
  dilution_levels <- c(10, 20, 40, 80, 160, 320, 640, 1280)
  aval_val_pool <- c(1:9, 91:99, 991:999, 1001:1009)

  # a small function to categorize aval values
  to_dilution_bucket <- function(x, levels = dilution_levels) {
    breaks <- c(0, levels, Inf)
    labels <- c(
      paste0("1:", levels),
      paste0(">", levels[length(levels)])
    )
    labels[findInterval(x, breaks)]
  }

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
      aval_val_pool,
      length(adatrept_idx),
      replace = TRUE
    )
  }

  # override AVAL for other titer PARAMCDs with proper doubling sequence values
  titer_idx <- which(main_tbl$PARAMCD %in% titer_paramcds & main_tbl$PARAMCD != "ADATREPT")
  if (length(titer_idx) > 0) {
    main_tbl$AVAL[titer_idx] <- sample(
      dilution_levels,
      length(titer_idx),
      replace = TRUE
    )
  }

  # row-wise small noise for numeric AVAL (non-titer only)
  noise <- stats::runif(n = nrow(main_tbl), min = -0.5, max = 0.5)

  main_tbl <- main_tbl |>
    dplyr::mutate(
      AVAL = dplyr::case_when(
        PARAMCD %in% titer_paramcds ~ AVAL, # keep titer AVAL as-is
        !is.na(ISSTRESN) & !PARAMCD %in% group_b_paramcd ~ round(ISSTRESN + noise, 2),
        .default = AVAL
      ),
      AVALC = dplyr::case_when(
        PARAMCD %in% titer_paramcds & !is.na(AVAL) ~ to_dilution_bucket(AVAL),
        !is.na(AVAL) & !PARAMCD %in% group_b_paramcd ~ as.character(floor(AVAL)),
        is.na(ISSTRESN) & !is.na(ISSTRESC) & !PARAMCD %in% group_b_paramcd ~ toupper(ISSTRESC),
        .default = AVALC
      ),
      AVALCAT1 = dplyr::case_when(
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL < 10 ~ "<10",
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL < 100 ~ "10 to <100",
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL < 1000 ~ "100 to <1000",
        PARAMCD == "ADATREPT" & !is.na(AVAL) & AVAL >= 1000 ~ ">=1000",
        !PARAMCD %in% group_b_paramcd &
          (!is.na(AVAL) & AVAL <= 0) |
          (!is.na(AVALC) & AVALC == "NEGATIVE") |
          (!is.na(AVALC) & grepl("^\\s*<", AVALC)) ~ "Negative / BLQ",
        !PARAMCD %in% group_b_paramcd & !is.na(AVAL) & AVAL > 0 & AVAL <= 2 ~ "Low Positive",
        !PARAMCD %in% group_b_paramcd & !is.na(AVAL) & AVAL > 2 ~ "High Positive",
        .default = AVALCAT1
      )
    )

  # NA fallback
  fallback_vals <- sample(aval_val_pool, nrow(main_tbl), replace = TRUE)
  main_tbl <- main_tbl |>
    dplyr::mutate(
      AVAL = dplyr::coalesce(AVAL, as.numeric(fallback_vals)),
      AVALC = dplyr::coalesce(AVALC, to_dilution_bucket(AVAL))
    )

  return(main_tbl)
}


# add IMEVFL column
add_adishum_col_funcs$imevfl <- function(main_tbl) {
  # Subject-level eligibility from IMFL
  subj <- main_tbl |>
    dplyr::group_by(USUBJID) |>
    dplyr::summarise(
      eligible = any(
        as.character(IMFL) == "Y",
        na.rm = TRUE
      ),
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
  visit_map <- add_adishum_col_funcs$sample_visits()

  main_tbl <- main_tbl |>
    dplyr::left_join(
      visit_map,
      by = "VISITNUM"
    ) |>
    dplyr::mutate(
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

# add ANL03FL, ANL04FL, ANL05FL, ANL06FL, ANL07FL,
#     ANL08FL, ANL09FL, ANL10FL, ANL11FL, ANL12FL columns
add_adishum_col_funcs$anl03fl_to_anl12fl <- function(main_tbl) {
  group_b_paramcd <- c("TMOSADAW", "ADADURW", "PSPDURW", "NABPOS", "NABNEG")

  pick_pct_03 <- 0.20
  pick_pct_07 <- 0.18

  # all treated subjects — one row per subject
  all_subj <- main_tbl |>
    dplyr::filter(!is.na(TRTSDT)) |>
    dplyr::distinct(USUBJID) |>
    dplyr::pull(USUBJID)

  # helper: pick a random fraction of a subject vector
  pick <- function(pool, pct) sample(pool, floor(length(pool) * pct))

  # ANL03FL: pick_pct_03 of all subjects; ANL07FL: pick_pct_07 of all subjects
  anl03_subj <- pick(all_subj, pick_pct_03)
  anl07_subj <- pick(all_subj, pick_pct_07)

  # derived flags — subsets of ANL03FL or ANL07FL subjects
  anl04_subj <- pick(anl03_subj, 0.35)
  anl05_subj <- pick(anl03_subj, 0.30)
  anl06_subj <- pick(anl03_subj, 0.25)
  anl08_subj <- pick(anl07_subj, 0.35)
  anl09_subj <- pick(anl07_subj, 0.30)
  anl10_subj <- pick(anl07_subj, 0.25)
  anl11_subj <- pick(anl03_subj, 0.40)
  anl12_subj <- pick(anl07_subj, 0.40)

  # join back: every row of a subject gets the same flag
  main_tbl <- main_tbl |>
    dplyr::mutate(
      ANL03FL = dplyr::if_else(USUBJID %in% anl03_subj, "Y", NA_character_),
      ANL04FL = dplyr::if_else(USUBJID %in% anl04_subj, "Y", NA_character_),
      ANL05FL = dplyr::if_else(USUBJID %in% anl05_subj, "Y", NA_character_),
      ANL06FL = dplyr::if_else(USUBJID %in% anl06_subj, "Y", NA_character_),
      ANL07FL = dplyr::if_else(USUBJID %in% anl07_subj, "Y", NA_character_),
      ANL08FL = dplyr::if_else(USUBJID %in% anl08_subj, "Y", NA_character_),
      ANL09FL = dplyr::if_else(USUBJID %in% anl09_subj, "Y", NA_character_),
      ANL10FL = dplyr::if_else(USUBJID %in% anl10_subj, "Y", NA_character_),
      ANL11FL = dplyr::if_else(USUBJID %in% anl11_subj, "Y", NA_character_),
      ANL12FL = dplyr::if_else(USUBJID %in% anl12_subj, "Y", NA_character_)
    )

  return(main_tbl)
}

# add atpt to my dataset
add_adishum_col_funcs$atpt_atptn <- function(main_tbl) {
  choices <- c("Pre-dose", "Dose", "8h Post-dose")
  is_summary <- tolower(main_tbl$PARCAT1) == "subject summary" |
    is.na(main_tbl$PARCAT1)

  atpt <- rep(NA_character_, nrow(main_tbl))
  atpt[!is_summary] <- sample(choices, sum(!is_summary, na.rm = TRUE), replace = TRUE)

  main_tbl$ATPT <- atpt
  main_tbl$ATPTN <- match(main_tbl$ATPT, choices) # Pre-dose=1, Dose=2, 8h Post-dose=3

  return(main_tbl)
}

# add INFUSPN, INFUSPRN, INFUSAN, INFUSARN,
#     INJECPN, INJECPRN, INJECAN, INJECARN columns
add_adishum_col_funcs$infus_injec_counts <- function(main_tbl) {
  # -- configurable ranges (change as needed) --
  infusTotalRange <- 1L:20L
  infusReactionRange <- 0L:5L
  injecTotalRange <- 1L:16L
  injecReactionRange <- 0L:4L

  subj <- main_tbl |>
    dplyr::distinct(USUBJID)

  n_subj <- nrow(subj)

  subj <- subj |>
    dplyr::mutate(
      INFUSPN = sample(infusTotalRange, n_subj, replace = TRUE),
      INFUSPRN = pmin(
        sample(infusReactionRange, n_subj, replace = TRUE),
        INFUSPN
      ),
      INFUSAN = sample(infusTotalRange, n_subj, replace = TRUE),
      INFUSARN = pmin(
        sample(infusReactionRange, n_subj, replace = TRUE),
        INFUSAN
      ),
      INJECPN = sample(injecTotalRange, n_subj, replace = TRUE),
      INJECPRN = pmin(
        sample(injecReactionRange, n_subj, replace = TRUE),
        INJECPN
      ),
      INJECAN = sample(injecTotalRange, n_subj, replace = TRUE),
      INJECARN = pmin(
        sample(injecReactionRange, n_subj, replace = TRUE),
        INJECAN
      )
    )

  main_tbl <- main_tbl |>
    dplyr::left_join(
      subj,
      by = "USUBJID"
    )

  return(main_tbl)
}

# add DTL (Drug Tolerance Limit) column - float in [0, 40] with noise to avoid perfect integers
add_adishum_col_funcs$dtl <- function(
  main_tbl,
  min_n = 0,
  max_n = 40
) {
  n <- nrow(main_tbl)
  main_tbl$DTL <- round(
    stats::runif(n, min = min_n, max = max_n) +
      stats::runif(n, min = -0.49, max = 0.49),
    digits = 2
  )
  main_tbl$DTL <- pmax(min_n, pmin(max_n, main_tbl$DTL))
  return(main_tbl)
}

# ensure required AVALC values are present per PARAMCD spec
add_adishum_col_funcs$pp_ensure_required_avalc <- function(
  main_tbl,
  min_rows = 10
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
    "ADAUND"
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

  # NABPOS/NABNEG: AVALC="Y" only for subjects who have ADATRE with AVALC="Y"
  adatre_y_subjects <- main_tbl |>
    dplyr::filter(PARAMCD == "ADATRE", AVALC == "Y") |>
    dplyr::distinct(USUBJID) |>
    dplyr::pull(USUBJID)

  main_tbl <- main_tbl |>
    dplyr::group_by(USUBJID, PARAMCD) |>
    dplyr::mutate(
      AVALC = dplyr::if_else(
        PARAMCD %in%
          c("NABPOS", "NABNEG") &
          USUBJID %in% adatre_y_subjects &
          dplyr::row_number() %in% sample(dplyr::n(), min(min_rows, dplyr::n())),
        "Y",
        AVALC
      )
    ) |>
    dplyr::ungroup()

  return(main_tbl)
}

# ensure required ANL09FL and ANAL05FL
add_adishum_col_funcs$pp_ensure_anl09fl_anl05fl <- function(main_tbl) {
  random_values <- 5

  eligible_subj <- main_tbl |>
    dplyr::filter(
      IMFL == "Y",
      PARAMCD == "ADATRE",
      AVALC == "Y"
    ) |>
    dplyr::distinct(USUBJID) |>
    dplyr::pull(USUBJID)

  pick05 <- sample(eligible_subj, min(random_values, length(eligible_subj)))
  pick09 <- sample(eligible_subj, min(random_values, length(eligible_subj)))

  main_tbl <- main_tbl |>
    dplyr::mutate(
      ANL05FL = dplyr::if_else(USUBJID %in% pick05, "Y", ANL05FL),
      ANL09FL = dplyr::if_else(USUBJID %in% pick09, "Y", ANL09FL)
    )

  return(main_tbl)
}

# ensure required values for ANL07FL and ANL10FL
add_adishum_col_funcs$pp_ensure_anl07fl_anl10fl <- function(main_tbl) {
  random_values <- 5

  # ensure ANL07FL has Y for some ADATRE positive subjects
  eligible_07 <- main_tbl |>
    dplyr::filter(
      SAFFL == "Y",
      PARAMCD == "ADATRE",
      AVALC == "Y"
    ) |>
    dplyr::distinct(USUBJID) |>
    dplyr::pull(USUBJID)

  pick07 <- sample(eligible_07, min(random_values, length(eligible_07)))

  main_tbl <- main_tbl |>
    dplyr::mutate(
      ANL07FL = dplyr::if_else(USUBJID %in% pick07, "Y", ANL07FL)
    )

  # ensure ANL10FL has Y — only from ANL07FL=Y subjects
  eligible_10 <- main_tbl |>
    dplyr::filter(
      SAFFL == "Y",
      PARAMCD == "ADATRE",
      AVALC == "Y",
      ANL07FL == "Y"
    ) |>
    dplyr::distinct(USUBJID) |>
    dplyr::pull(USUBJID)

  pick10 <- sample(eligible_10, min(random_values, length(eligible_10)))

  main_tbl <- main_tbl |>
    dplyr::mutate(
      ANL10FL = dplyr::if_else(USUBJID %in% pick10, "Y", ANL10FL)
    )

  return(main_tbl)
}

# ensure that for ANL02FL only the last row per (USUBJID, AVISIT, ATPT) group is "Y"
add_adishum_col_funcs$pp_ensure_anl02fl <- function(main_tbl) {
  main_tbl <- main_tbl |>
    dplyr::arrange(USUBJID, AVISIT, ATPT) |>
    dplyr::group_by(USUBJID, AVISIT, ATPT) |>
    dplyr::mutate(
      ANL02FL = dplyr::if_else(
        dplyr::row_number() == dplyr::n(),
        "Y",
        NA_character_
      )
    ) |>
    dplyr::ungroup()

  return(main_tbl)
}

# ensure SUADAST rows cover all 9 AVISIT+ATPT combos with min POSITIVE and NEGATIVE per combo
add_adishum_col_funcs$pp_ensure_suadast_avalc <- function(
  main_tbl,
  n_day1 = 40,
  n_day3 = 40,
  n_per_atpt = 12,
  pos_range = 2:6,
  neg_range = 1:3
) {
  main_tbl$AVALC <- as.character(main_tbl$AVALC)
  atpt_choices <- c("Pre-dose", "Dose", "8h Post-dose")

  suadast_idx <- which(main_tbl$PARAMCD == "SUADAST")

  # step 1: pick ~n_day1 existing SUADAST rows, update AVISIT to Day 1; ~n_day3 to Day 3
  day1_idx <- sample(suadast_idx, min(n_day1, length(suadast_idx)))
  remaining <- setdiff(suadast_idx, day1_idx)
  day3_idx <- sample(remaining, min(n_day3, length(remaining)))

  main_tbl$AVISIT[day1_idx] <- "Day 1"
  main_tbl$AVISITN[day1_idx] <- 1L
  main_tbl$AVISIT[day3_idx] <- "Day 3"
  main_tbl$AVISITN[day3_idx] <- 3L

  # step 2: assign ATPT within day1 and day3 rows so each choice gets ~n_per_atpt rows
  for (idx in list(day1_idx, day3_idx)) {
    n <- length(idx)
    atpt_vec <- sample(rep(atpt_choices, each = n_per_atpt)[seq_len(n)])
    main_tbl$ATPT[idx] <- atpt_vec
    main_tbl$ATPTN[idx] <- match(atpt_vec, atpt_choices)
  }

  # step 3: for each non-NA AVISIT+ATPT combo, insert min 2 to max 6 POSITIVE and min 1 to max 3 NEGATIVE
  valid_idx <- which(
    main_tbl$PARAMCD == "SUADAST" &
      !is.na(main_tbl$AVISIT) &
      !is.na(main_tbl$ATPT)
  )

  combos <- unique(main_tbl[valid_idx, c("AVISIT", "ATPT")])

  for (i in seq_len(nrow(combos))) {
    cidx <- valid_idx[
      main_tbl$AVISIT[valid_idx] == combos$AVISIT[i] &
        main_tbl$ATPT[valid_idx] == combos$ATPT[i]
    ]
    n_pos <- min(sample(pos_range, 1), length(cidx))
    n_neg <- min(sample(neg_range, 1), length(cidx) - n_pos)
    chosen <- sample(cidx, n_pos + n_neg)
    main_tbl$AVALC[chosen[seq_len(n_pos)]] <- "POSITIVE"
    main_tbl$AVALC[chosen[seq(n_pos + 1, n_pos + n_neg)]] <- "NEGATIVE"
  }

  return(main_tbl)
}

# add tests to my dataset
add_adishum_col_funcs$dataset_tests <- function(main_tbl) {
  rule1 <- main_tbl |>
    count(
      USUBJID,
      VISITNUM,
      PARAMCD
    ) |>
    filter(
      n > 1
    )

  if (nrow(rule1) != 0) {
    message(paste(rep("-", 80)))
    message("FALSE -- Each USUBJID should have 1 row per VISITNUM per PARAMCD")
    message(paste(rep("-", 80)))
  }

  rule2 <- main_tbl |>
    dplyr::select(
      dplyr::starts_with("ANL")
    ) |>
    lapply(\(x) any(x == "Y", na.rm = TRUE))

  if (!all(unlist(rule2))) {
    message(paste(rep("-", 80)))
    message('FALSE -- No "Y" value in one or more ANL__FL columns')
    message(paste(rep("-", 80)))
  }

  rule3 <- main_tbl |>
    filter(
      PARAMCD == "ADATREPT"
    ) |>
    count(
      AVALCAT1
    )

  avalcat_values <- c("<10", "10 to <100", "100 to <1000", ">=1000")

  if (
    !all(rule3$AVALCAT1 %in% avalcat_values) &&
      any(rule3$n != 0)
  ) {
    message(paste(rep("-", 80)))
    message(r"{FALSE -- AVALCAT1 has some values when PARAMCD == "ADATREPT"}")
    message(paste(rep("-", 80)))
  }

  rule4 <- main_tbl |>
    dplyr::distinct(TRT01A, TRT01AN)

  if (any(duplicated(rule4$TRT01A))) {
    message(paste(rep("-", 80)))
    message("FALSE -- Treatment columns have different values")
    message(paste(rep("-", 80)))
  }

  anl_flags <- paste0("ANL", sprintf("%02d", 3:12), "FL")

  rule5 <- main_tbl |>
    dplyr::count(
      USUBJID,
      dplyr::across(dplyr::all_of(anl_flags))
    ) |>
    dplyr::count(USUBJID) |>
    dplyr::filter(n > 1)

  if (nrow(rule5) != 0) {
    message(paste(rep("-", 80)))
    message("FALSE -- ANL03FL to ANL12FL are not consistent at subject level")
    message(paste(rep("-", 80)))
  }

  rule6 <- main_tbl |>
    dplyr::filter(TRT01A == "Placebo") |>
    dplyr::filter(IMFL == "Y")

  if (nrow(rule6) != 0) {
    message(paste(rep("-", 80)))
    message("FALSE -- Placebo subjects should have IMFL = N")
    message(paste(rep("-", 80)))
  }

  return(TRUE)
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
    "ISBDAGNT",
    "ISSTRESC",
    "ISSTRESN",
    "ISNAM",
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

  adsl <- pharmaverseadamjnj::adsl

  gen <- raw |>
    dplyr::select(
      dplyr::all_of(is_ada_vars)
    ) |>
    add_adishum_col_funcs$sync_adsl_vars(
      adsl_vars,
      adsl
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
    add_adishum_col_funcs$anl03fl_to_anl12fl()

  # ensure required ATPT column is added
  gen <- gen |>
    add_adishum_col_funcs$atpt_atptn()

  # add INFUSPN, INFUSPRN, INFUSAN, INFUSARN,
  #     INJECPN, INJECPRN, INJECAN, INJECARN columns - refer function
  gen <- gen |>
    add_adishum_col_funcs$infus_injec_counts()

  # add DTL (Drug Tolerance Limit) column - refer function
  gen <- gen |>
    add_adishum_col_funcs$dtl()

  # POST PROCESSING for inserting random values to meet standards -------------
  # ensure required AVALC values are present per PARAMCD spec
  gen <- gen |>
    add_adishum_col_funcs$pp_ensure_required_avalc()

  # ensure required ANL09FL and ANAL05FL are there
  gen <- gen |>
    add_adishum_col_funcs$pp_ensure_anl09fl_anl05fl()

  # ensure required ANL09FL and ANAL05FL are there
  gen <- gen |>
    add_adishum_col_funcs$pp_ensure_anl07fl_anl10fl()

  # ensure that only last row per (USUBJID, AVISIT, ATPT) gets "Y" for ANL02FL
  gen <- gen |>
    add_adishum_col_funcs$pp_ensure_anl02fl()

  # ensure that we have a few rows where AVALC="Y" when PARAMCD="SUADAST"
  gen <- gen |>
    add_adishum_col_funcs$pp_ensure_suadast_avalc()

  # sync adsl vars at the end to fix any corruption from row sampling
  gen <- gen |>
    add_adishum_col_funcs$sync_adsl_vars(adsl_vars, adsl)

  # run basic tests: ----------------------------------
  # ensure Dataset is valid for working
  gen |>
    add_adishum_col_funcs$dataset_tests()

  # Handle NA values and convert characters to factors
  gen <- gen |>
    df_na(char_as_factor = TRUE)

  # Define additional labels for new variables
  additional_labels <- list(
    IMFL = "Immunogenicity Population Flag",
    PARCAT1 = "Parameter Category 1",
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
    ANL11FL = "Analysis Flag 11-grade >=3 Inf SR",
    ANL12FL = "Analysis Flag 12-grade >=3 Inj SR",
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
    TRTA = "Actual Treatment",
    TRTSDT = "Date of First Exposure to Treatment",
    TRTSDTM = "Datetime of First Exposure to Treatment",
    TRT01P = "Planned Treatment for Period 01",
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01A = "Actual Treatment for Period 01",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    ATPT = "Analysis Timepoint",
    ATPTN = "Analysis Timepoint (N)",
    INFUSPN = "Number of Placebo infusions",
    INFUSPRN = "Placebo infusions with reaction",
    INFUSAN = "Number of Active Drug infusions",
    INFUSARN = "Active Drug infusions with reaction",
    INJECPN = "Number of Placebo injections",
    INJECPRN = "Placebo injections with reaction",
    INJECAN = "Number of Active Drug injections",
    INJECARN = "Active Drug injections with reaction",
    DTL = "Drug Tolerance Limit"
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
