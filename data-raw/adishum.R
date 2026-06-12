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

# function for adding PARQUAL col
add_parqual_col <- function(main_tbl) {
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
add_param_n_paramcd_col <- function(main_tbl) {
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
add_aval_avalc_avalcat1_col <- function(main_tbl) {
  # add empty columns
  main_tbl <- main_tbl |>
    mutate(
      AVAL  = as.numeric(NA),
      AVALC = as.character(NA),
      AVALCAT1 = as.character(NA)
    )

  # classification using stringi
  is_titer  <- stringi::stri_endswith_fixed(str = main_tbl$PARAMCD, pattern = "T")
  is_nabposneg <- main_tbl$PARAMCD %in% c("NABPOS", "NABNEG")
  is_result_text <- stringi::stri_detect_regex(
    str = main_tbl$PARAM,
    pattern = "result|positive|negative|undetermined",
    opts_regex = stringi::stri_opts_regex(case_insensitive = TRUE)
  )

  # ensure titers take priority
  is_result <- ifelse(
    is_titer,
    FALSE,
    (is_nabposneg | is_result_text)
  )

  # simple titer pool
  titer_values <- c(10, 20, 40, 80, 160, 320)

  # populate AVALC for titers (include some left-censored "<" values)
  n_t <- sum(is_titer)
  if (n_t > 0) {
    vals <- sample(titer_values, size = n_t, replace = TRUE)
    censored <- runif(n_t) < 0.20   # 20% left-censored example
    main_tbl$AVALC[which(is_titer)] <- paste0(ifelse(censored, "<", ""), vals)
  }

  # populate AVALC for result-like rows
  n_r <- sum(is_result)
  if (n_r > 0) {
    main_tbl$AVALC[which(is_result)] <- sample(
      c("Positive", "Negative", "Undetermined"),
      size = n_r,
      replace = TRUE,
      prob = c(0.25, 0.70, 0.05)
    )
  }

  # extract numeric part and left-censor marker using stringi
  numpart <- stringi::stri_extract_first_regex(str = main_tbl$AVALC, pattern = "[0-9]+")
  numval  <- as.numeric(numpart)
  is_lt   <- stringi::stri_startswith_fixed(str = main_tbl$AVALC, pattern = "<")
  is_digits <- stringi::stri_detect_regex(str = main_tbl$AVALC, pattern = "^[0-9]+$")
  
  # derive numeric AVAL and round to 2 significant digits
  main_tbl <- main_tbl |>
    mutate(
      AVAL = dplyr::case_when(
        !is.na(AVALC) & is_lt ~ numval / 2,
        !is.na(AVALC) & is_digits ~ numval,
        !is.na(AVALC) & stringi::stri_trans_tolower(AVALC) == "positive" ~ 1,
        !is.na(AVALC) & stringi::stri_trans_tolower(AVALC) == "negative" ~ 0,
        TRUE ~ as.numeric(NA)
      ),
      AVAL = ifelse(!is.na(AVAL), signif(AVAL, 2), NA_real_)
    )

  # simple categorical bins for titers and passthrough for results
  main_tbl <- main_tbl |>
    mutate(
      AVALCAT1 = dplyr::case_when(
        !is.na(AVAL) & is_titer & AVAL <= 20  ~ "Negative",
        !is.na(AVAL) & is_titer & AVAL <= 80  ~ "Low",
        !is.na(AVAL) & is_titer & AVAL <= 320 ~ "Medium",
        !is.na(AVAL) & is_titer & AVAL > 320  ~ "High",
        !is.na(AVALC) & is_result ~ AVALC,  # use original text for result rows
        TRUE ~ NA_character_
      )
    )
  
  return(main_tbl)
}

# core function ------------------------------
# Generate adishum dataset
gen_adishum <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # get source data
  raw <- pharmaverseadam::adis_vaccine

  gen <- raw |> 
    select(
      STUDYID,
      USUBJID
    )

  sdtm_tbl <- pharmaversesdtm::is_ada |> 
    left_join(
      pharmaverseadamjnj::adsl,
      by = "USUBJID"
    ) |> 
    select(
      IMFL, 
      SAFFL,
      ISBDAGNT,
      TRTSDT, 
      TRTSDTM,
      TRT01P, 
      TRT01PN, 
      TRT01A, 
      TRT01AN, 
      AGE, 
      SEX, 
      RACE, 
      SITEID, 
      SUBJID,
      ISTESTCD
    )

  # add PARQUAL column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_parqual_col()

  # add PARAM and PARAMCD column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_param_n_paramcd_col()

  # add AVAL, AVALC and AVALCAT1 column - refer function
  sdtm_tbl <- sdtm_tbl |> 
    add_aval_avalc_avalcat1_col()




}


# run function ------------------------------
# generate dataset
adishum <- gen_adishum()
