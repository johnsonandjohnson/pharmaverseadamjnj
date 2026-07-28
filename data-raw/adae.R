# Generate ADAE dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(formatters)
library(forcats)
library(admiral)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADAE dataset
gen_adae <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  # Get source data
  raw <- pharmaverseadam::adae

  # Source ADSL to get treatment duration information
  source(file.path("data-raw", "adsl.R"))

  # Source ADISHUM to get ADATRES
  source(file.path("data-raw", "adishum.R"))

  # Create a mapping of USUBJID to their appropriate ACAT1 category based on TRTEDY
  subject_acat1_map <- adsl |>
    dplyr::mutate(
      ACAT1 = dplyr::case_when(
        TRTEDY <= 90 ~ "Within 3 months",
        TRTEDY > 90 & TRTEDY <= 180 ~ "4 to 6 months",
        TRTEDY > 180 & TRTEDY <= 270 ~ "7 to 9 months",
        TRTEDY > 270 & TRTEDY <= 365 ~ "10 to 12 months",
        TRTEDY > 365 ~ "Beyond 13 months"
      )
    ) |>
    dplyr::select(USUBJID, ACAT1)

  gen <- raw

  # Update specific records as per requirements:

  # a.  First records Placebo and Treatment for Female Narrow
  gen <- gen |>
    mutate(
      AETERM = ifelse(USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1, "ABNORMAL UTERINE BLEEDING", AETERM),
      AELLT = ifelse(
        USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1,
        "DYSFUNCTIONAL UTERINE BLEEDING",
        AELLT
      ),
      AEDECOD = ifelse(USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1, "ABNORMAL UTERINE BLEEDING", AEDECOD),
      AEHLT = ifelse(
        USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1,
        "MENSTRUATION AND UTERINE BLEEDING NEC",
        AEHLT
      ),
      AEHLGT = ifelse(
        USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1,
        "MENSTRUAL CYCLE AND UTERINE BLEEDING DISORDERS",
        AEHLGT
      ),
      AEBODSYS = ifelse(
        USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1,
        "REPRODUCTIVE SYSTEM AND BREAST DISORDERS",
        AEBODSYS
      ),
      AESOC = ifelse(
        USUBJID %in% c("01-701-1015", "01-701-1034") & AESEQ == 1,
        "REPRODUCTIVE SYSTEM AND BREAST DISORDERS",
        AESOC
      )
    )

  # b.  First records Placebo and Treatment for Male Narrow
  gen <- gen |>
    mutate(
      AETERM = ifelse(USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1, "ERECTILE DYSFUNCTION", AETERM),
      AELLT = ifelse(USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1, "ERECTILE DISTURBANCE", AELLT),
      AEDECOD = ifelse(USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1, "ERECTILE DYSFUNCTION", AEDECOD),
      AEHLT = ifelse(
        USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1,
        "ERECTION AND EJACULATION CONDITIONS AND DISORDERS",
        AEHLT
      ),
      AEHLGT = ifelse(
        USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1,
        "SEXUAL FUNCTION AND FERTILITY DISORDERS",
        AEHLGT
      ),
      AEBODSYS = ifelse(
        USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1,
        "REPRODUCTIVE SYSTEM AND BREAST DISORDERS",
        AEBODSYS
      ),
      AESOC = ifelse(
        USUBJID %in% c("01-701-1023", "01-701-1028") & AESEQ == 1,
        "REPRODUCTIVE SYSTEM AND BREAST DISORDERS",
        AESOC
      )
    )

  # c.  First records Placebo and Treatment for Female Broad
  gen <- gen |>
    mutate(
      AETERM = ifelse(USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1, "BLEEDING ANOVULATORY", AETERM),
      AELLT = ifelse(
        USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1,
        "ANOVULAR DYSFUNCTIONAL UTERINE BLEEDING",
        AELLT
      ),
      AEDECOD = ifelse(USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1, "BLEEDING ANOVULATORY", AEDECOD),
      AEHLT = ifelse(
        USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1,
        "FEMALE GONADAL FUNCTION DISORDERS",
        AEHLT
      ),
      AEHLGT = ifelse(
        USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1,
        "ENDOCRINE DISORDERS OF GONADAL FUNCTION",
        AEHLGT
      ),
      AEBODSYS = ifelse(USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1, "ENDOCRINE DISORDERS", AEBODSYS),
      AESOC = ifelse(USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1, "ENDOCRINE DISORDERS", AESOC),
      TRTEMFL = ifelse(USUBJID %in% c("01-701-1363", "01-701-1111") & AESEQ == 1, "Y", TRTEMFL)
    )

  # d.  First records Placebo and Treatment for Male Broad
  gen <- gen |>
    mutate(
      AETERM = ifelse(
        USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1,
        "DISTURBANCE IN SEXUAL AROUSAL",
        AETERM
      ),
      AELLT = ifelse(USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1, "SEXUAL AROUSAL DECREASED", AELLT),
      AEDECOD = ifelse(
        USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1,
        "DISTURBANCE IN SEXUAL AROUSAL",
        AEDECOD
      ),
      AEHLT = ifelse(USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1, "SEXUAL AROUSAL DISORDERS", AEHLT),
      AEHLGT = ifelse(
        USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1,
        "SEXUAL DYSFUNCTIONS, DISTURBANCES AND GENDER IDENTITY DISORDERS",
        AEHLGT
      ),
      AEBODSYS = ifelse(USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1, "PSYCHIATRIC DISORDERS", AEBODSYS),
      AESOC = ifelse(USUBJID %in% c("01-701-1392", "01-701-1097") & AESEQ == 1, "PSYCHIATRIC DISORDERS", AESOC)
    )

  gen <- dplyr::mutate(
    gen,
    # Add times to ASTDTM
    ASTDTM = {
      # Extract the date part from ASTDTM
      date_parts <- sub(" UTC", "", ASTDTM)
      # Generate random hours (0-23) for each row
      # Format time strings (00:00:00)
      times <- sprintf("00:00:00")
      # Combine date with random time and UTC timezone
      paste(date_parts, times, "UTC")
    },
    AETOXGR = as.factor(sample(seq(0, 5), dplyr::n(), replace = TRUE)),
    AETOXGRN = as.numeric(as.character(AETOXGR)),
    AEACN = as.factor(sample(
      c(
        "DOSE NOT CHANGED",
        "NOT APPLICABLE",
        "DRUG WITHDRAWN",
        "DOSE REDUCED",
        "DRUG INTERRUPTED",
        "DOSE RATE REDUCED",
        "DOSE INCREASED",
        "UNKNOWN"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AEACNS1 = as.factor(sample(
      c(
        "DOSE NOT CHANGED",
        "NOT APPLICABLE",
        "DRUG WITHDRAWN",
        "DOSE REDUCED",
        "DRUG INTERRUPTED",
        "DOSE RATE REDUCED",
        "DOSE INCREASED",
        "UNKNOWN"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AEACNS2 = as.factor(sample(
      c(
        "DOSE NOT CHANGED",
        "NOT APPLICABLE",
        "DRUG WITHDRAWN",
        "DOSE REDUCED",
        "DRUG INTERRUPTED",
        "DOSE RATE REDUCED",
        "DOSE INCREASED",
        "UNKNOWN"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AESEV = dplyr::case_when(
      AESEV == "MILD" ~ "Mild",
      AESEV == "MODERATE" ~ "Moderate",
      AESEV == "SEVERE" ~ "Severe"
    ),
    DOSEDY = as.numeric(37),
    DOSS1DY = as.numeric(38),
    DOSS2DY = as.numeric(39),
    DOSEU = as.factor("mg"),
    DOSS1U = as.factor("mcg"),
    DOSS2U = as.factor("mL"),
    DOSEON = as.numeric(10),
    DOSS1ON = as.numeric(11),
    DOSS2ON = as.numeric(12),
    AECONTRT = as.factor(sample(c("N", "Y", "U"), dplyr::n(), replace = TRUE)),
    CQ01NAM = as.factor(sample(c("Seizure", NA), dplyr::n(), replace = TRUE)),
    CQ02NAM = as.factor(sample(c("Skin rash", NA), dplyr::n(), replace = TRUE)),
    CQ03NAM = as.factor(sample(
      c("Hypothyroidism", NA),
      dplyr::n(),
      replace = TRUE
    )),
    AESMIE = "Y",
    AESER = as.factor(sample(c("N", "Y"), dplyr::n(), replace = TRUE)),
    AEREL = as.factor(dplyr::case_when(
      AEREL == "PROBABLE" ~ "RELATED",
      AEREL == "REMOTE" ~ "RELATED",
      AEREL == "POSSIBLE" ~ "RELATED",
      AEREL == "NONE" ~ "NOT RELATED",
      is.na(AEREL) ~ NA_character_
    )),
    AERELS1 = as.factor(sample(
      c(
        "RELATED",
        "NOT RELATED"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    AERELS2 = as.factor(sample(
      c(
        "RELATED",
        "NOT RELATED"
      ),
      dplyr::n(),
      replace = TRUE
    )),
    RELGR1 = as.factor(dplyr::case_when(
      if_any(starts_with("AERELS"), ~ . == "RELATED") ~ "RELATED",
      if_all(starts_with("AERELS"), ~ !is.na(.)) ~ "NOT RELATED",
      TRUE ~ NA_character_
    )),
    AEDRGS1 = "Study agent 1",
    AEDRGS2 = "Study agent 2",
    AEDECOD = forcats::fct_relabel(AEDECOD, stringr::str_to_sentence), # Convert AEDECOD levels to sentence
    AEBODSYS = forcats::fct_relabel(AEBODSYS, stringr::str_to_sentence) # Convert AEBODSYS levels to sentence
  )

  # Join with subject_acat1_map to assign ACAT1 based on subject's treatment duration
  gen <- dplyr::left_join(gen, subject_acat1_map, by = "USUBJID")
  gen$ACAT1 <- as.factor(gen$ACAT1)

  # Apply derivations
  gen <- gen |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID),
        order = exprs(STUDYID, !is.na(ASTDT), ASTDT, AESEQ),
        new_var = AOCCFL,
        mode = "first"
      ),
      filter = TRTEMFL == "Y"
    ) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, AEDECOD),
        order = exprs(STUDYID, !is.na(ASTDT), ASTDT, AESEQ),
        new_var = AOCCPFL,
        mode = "first"
      ),
      filter = TRTEMFL == "Y"
    ) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, AEBODSYS),
        order = exprs(STUDYID, !is.na(ASTDT), ASTDT, AESEQ),
        new_var = AOCCSFL,
        mode = "first"
      ),
      filter = TRTEMFL == "Y"
    )

  te_flags <- gen |>
    dplyr::filter(TRTEMFL == "Y") |>
    derive_var_extreme_flag(
      new_var = AOCTIFL,
      by_vars = exprs(USUBJID),
      order = exprs(dplyr::desc(AETOXGRN), ASTDY, AESEQ),
      mode = "first",
      false_value = "N"
    ) |>
    dplyr::select(USUBJID, AESEQ, AOCTIFL)

  gen <- gen |>
    dplyr::left_join(te_flags, by = c("USUBJID", "AESEQ"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), c("USUBJID", "TRTEDY"))

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c(
    "TRT01A",
    "SAFFL",
    "AGE",
    "SEX",
    "RACE",
    "STUDYID",
    "AGEGR1",
    "TRTEDY",
    "TRT01P"
  )

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl |>
    select(USUBJID, all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    message("Dropping shared vars from raw: ", paste(shared, collapse = ", "))
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID")
  gen <- mutate(
    gen,
    months = (TRTEDY + 30) / 30.4375,
    ACAT1 = case_when(
      months <= 3 ~ "Within 3 months",
      months > 3 & months <= 6 ~ "4 to 6 months",
      months > 6 & months <= 9 ~ "7 to 9 months",
      months > 9 & months <= 12 ~ "10 to 12 months",
      months > 12 ~ "Beyond 13 months",
      .default = NA_character_
    ),
    ACAT1 = factor(
      ACAT1,
      levels = c(
        "Within 3 months",
        "4 to 6 months",
        "7 to 9 months",
        "10 to 12 months",
        "Beyond 13 months"
      )
    )
  )
  gen <- select(gen, -months)

  # Add TRDISCFL variable: "Y" if AEACN = "DRUG WITHDRAWN", null otherwise
  gen <- gen |>
    mutate(
      TRDISCFL = case_when(
        AEACN == "DRUG WITHDRAWN" ~ "Y",
        AEACN == "MULTIPLE" & if_any(starts_with("AEACNS"), ~ . == "DRUG WITHDRAWN") ~ "Y",
        TRUE ~ NA_character_
      )
    )

  gen <- gen |>
    mutate(
      AESHOSPP = ifelse(AESHOSP == "Y", "Y", NA_character_),
      AESHOSPR = ifelse(AESHOSP == "Y", "Y", NA_character_)
    )

  gen <- gen |>
    mutate(
      AESCAT = case_when(
        AESOC == "GENERAL DISORDERS AND ADMINISTRATION SITE CONDITIONS" ~ "INFUSION RELATED REACTION",
        .default = "NONE OF THE ABOVE"
      ),
      AESCAT = forcats::fct_relabel(AESCAT, stringr::str_to_sentence), # Convert AESCAT levels to sentence
    )

  # Derive CQ01/02/03 and SMQ01/02/03 names per mapping from AEDECOD
  gen <- gen |>
    mutate(
      .AEDECOD_UP = toupper(AEDECOD),
      CQ01NAM = case_when(
        .AEDECOD_UP == "HYPERTENSION" ~ "Hypertension",
        .default = NA_character_
      ),
      SMQ01NAM = case_when(
        .AEDECOD_UP == "HYPERTENSION" ~ "Hypertension",
        .default = NA_character_
      ),
      CQ02NAM = case_when(
        .AEDECOD_UP %in% c("PRURITUS", "ERYTHEMA", "RASH") ~ "Sensitivity",
        .default = NA_character_
      ),
      SMQ02NAM = case_when(
        .AEDECOD_UP %in% c("PRURITUS", "ERYTHEMA", "RASH") ~ "Hypersensitivity",
        .default = NA_character_
      ),
      CQ03NAM = case_when(
        .AEDECOD_UP == "DIZZINESS" ~ "Hearing disorders",
        .default = NA_character_
      ),
      SMQ03NAM = case_when(
        .AEDECOD_UP == "DIZZINESS" ~ "Hearing and vestibular disorders",
        .default = NA_character_
      )
    ) |>
    select(-.AEDECOD_UP)

  # Derive AOCTxxFL variables
  gen <- gen |>
    arrange(USUBJID, desc(AETOXGRN), ASTDT) |>
    group_by(USUBJID) |>
    mutate(
      across(
        starts_with("CQ"),
        ~ {
          target_row <- match(
            TRUE,
            TRTEMFL == "Y" & !is.na(.x)
          )
          if_else(
            row_number() == target_row,
            "Y",
            NA_character_
          )
        },
        .names = "{gsub('CQ([0-9]+)NAM', 'AOCT\\\\1FL', .col)}"
      )
    ) |>
    ungroup()

  # Derive AOCSxxFL variables
  gen <- gen |>
    arrange(USUBJID, desc(ASEVN), ASTDT) |>
    group_by(USUBJID) |>
    mutate(
      across(
        starts_with("CQ"),
        ~ {
          target_row <- match(
            TRUE,
            TRTEMFL == "Y" & !is.na(.x)
          )
          if_else(
            row_number() == target_row,
            "Y",
            NA_character_
          )
        },
        .names = "{gsub('CQ([0-9]+)NAM', 'AOCS\\\\1FL', .col)}"
      )
    ) |>
    ungroup()

  # Add labels
  additional_labels <- list(
    SAFFL = "Safety Population Flag",
    AESER = "Serious Event",
    ACAT1 = "Analysis Category 1",
    AETOXGR = "Standard Toxicity Grade",
    AETOXGRN = "Standard Toxicity Grade (N)",
    AOCTIFL = "1st TE Max Toxicity Grade Flag",
    DOSEDY = "Day of Study Drug",
    DOSEU = "Treatment Dose Units",
    DOSEON = "Treatment Dose at Record Start",
    AOCCFL = "1st Occurrence within Subject Flag",
    AOCCPFL = "1st Occurrence within Pref Term Flag",
    AOCCSFL = "1st Occurrence of SOC Flag",
    CQ01NAM = "Customized Query 01 Name",
    CQ02NAM = "Customized Query 02 Name",
    CQ03NAM = "Customized Query 03 Name",
    SMQ01NAM = "Standardized MedDRA Query 01 Name",
    SMQ02NAM = "Standardized MedDRA Query 02 Name",
    SMQ03NAM = "Standardized MedDRA Query 03 Name",
    AECONTRT = "Concomitant or Additional Trtmnt Given",
    AESMIE = "Other Medically Important Serious Event",
    TRDISCFL = "Treatment Discontinued Flag",
    AESHOSPP = "Prolongs Hospitalization",
    AESHOSPR = "Requires Hospitalization",
    AOCT01FL = "1st AESI Max Tox. Grade 01 Occur Flag",
    AOCT02FL = "1st AESI Max Tox. Grade 02 Occur Flag",
    AOCT03FL = "1st AESI Max Tox. Grade 03 Occur Flag",
    AOCS01FL = "1st AESI Max Sev./Int. 01 Occur. Flag",
    AOCS02FL = "1st AESI Max Sev./Int. 02 Occur. Flag",
    AOCS03FL = "1st AESI Max Sev./Int. 03 Occur. Flag",
    AESCAT = "Adverse Event Category",
    AERELS1 = "Causality - Sponsor Study Treatment 1",
    AERELS2 = "Causality - Sponsor Study Treatment 2",
    AEACNS1 = "Action Taken - Sponsor Study Treatment 1",
    AEACNS2 = "Action Taken - Sponsor Study Treatment 2",
    RELGR1 = "Pooled Causality Group 1",
    AEDRGS1 = "Sponsor Study Treatment 1",
    AEDRGS2 = "Sponsor Study Treatment 2",
    DOSS1DY = "Day of Study Drug of study Agent 1",
    DOSS2DY = "Day of Study Drug of study Agent 2",
    DOSS1U = "Trt Dose Units for study Agent 1",
    DOSS2U = "Trt Dose Units for study Agent 2",
    DOSS1ON = "Treatment Dose for study Agent 1",
    DOSS2ON = "Treatment Dose for study Agent 2",
    ADATRES = "Treatment-emergent ADA Subject Status"
  )

  # Join ADATRES from ADISHUM - one row per subject, POSITIVE takes priority
  adatres_map <- adishum |>
    dplyr::filter(!is.na(ADATRES)) |>
    dplyr::mutate(priority_flag = dplyr::if_else(ADATRES == "POSITIVE", 1L, 2L)) |>
    dplyr::arrange(USUBJID, priority_flag) |>
    dplyr::distinct(USUBJID, .keep_all = TRUE) |>
    dplyr::select(USUBJID, ADATRES)

  gen <- gen |>
    dplyr::left_join(adatres_map, by = "USUBJID")

  # Arrange final data
  gen <- gen |> arrange(USUBJID, AESEQ)

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

# Generate the dataset
adae <- gen_adae()
