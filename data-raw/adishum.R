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


    





}


# run function ------------------------------
# generate dataset
adishum <- gen_adishum()
