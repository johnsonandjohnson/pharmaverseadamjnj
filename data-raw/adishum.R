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
      SUBJID
    )
  
  sdtm_tbl <- sdtm_tbl |> 
    mutate(
      PARQUAL = if_else(
        is.na(ISBDAGNT), 
        NA_character_,
        ISBDAGNT |> 
          str_trim() |> 
          str_to_upper()
    )
  )

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
    "ADATREPT", "Treatment-Enhanced ADA, Titer (duplicate entry)",          "ADA_BAB",
    "ADANTRE",  "Treatment-Enhanced Neutralizing ADA, Result",              "ADA_NAB",
    "ADAPSP",   "ADA Persistent Positive, Result",                          "ADA_BAB",
    "ADATSP",   "Neutralizing ADA Persistent Positive",                     "ADA_NAB",
    "ADAUND",   "ADA Undetermined, Result",                                 "ADA_BAB",
    "TMOSADAW", "Time to Onset of ADA (weeks)",                             NA_character_,
    "ADADURW",  "Duration of ADA Response (weeks)",                         NA_character_,
    "PSPDURW",  "Duration of Persistent Positive ADA (weeks)",              NA_character_,
    "NABPOS",   "Neutralizing ADA Positive",                                "ADA_NAB",
    "NABNEG",   "Neutralizing ADA Negative",                                "ADA_NAB"
  )

  random_sample_paramcd <- sample(1:20, nrow(sdtm_tbl), replace = TRUE)

  temp_paramcd <- temp_paramcd[
    random_sample_paramcd,
  ]

  sdtm_tbl <- sdtm_tbl |> 
    cbind(temp_paramcd)



    





}


# run function ------------------------------
# generate dataset
adishum <- gen_adishum()
