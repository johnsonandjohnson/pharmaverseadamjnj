# Generate ADEX dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(formatters)
library(forcats)
library(stringr)

# Source utility functions
source(file.path("data-raw", "helpers.R"))
# Ensure ADAE is available for cross-domain derivations
source(file.path("data-raw", "adae.R"))

# Source ADISHUM to get ADATRES
source(file.path("data-raw", "adishum.R"))

# Generate ADEX dataset
gen_adex <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  # Get source data
  raw <- pharmaverseadam::adex

  raw <- dplyr::filter(raw, PARAMCD == "DOSE")

  gen <- raw

  gen$TRT01P <- as.factor(gen$TRT01P)
  gen$TRT01A <- as.factor(gen$TRT01A)

  gen <- dplyr::mutate(
    gen,
    ECDOSE = dplyr::case_when(
      EXDOSE == 54 ~ 7.5,
      EXDOSE == 81 ~ 15,
      .default = EXDOSE
    ),
    ECOCCUR = factor(
      dplyr::case_when(
        !is.na(ECDOSE) ~ "Y",
        .default = "N"
      ),
      levels = c("Y", "N")
    ),
    EXDOSE = dplyr::case_when(
      EXDOSE == 54 ~ 7.5,
      EXDOSE == 81 ~ 15,
      .default = EXDOSE
    ),
    ATRT = factor(
      dplyr::case_when(
        TRT01P == "Xanomeline High Dose" ~ "XANOMELINE",
        TRT01P == "Xanomeline Low Dose" ~ "XANOMELINE",
        TRT01P == "Placebo" ~ "PLACEBO"
      ),
      levels = c(
        "XANOMELINE",
        "PLACEBO"
      )
    ),
    DAEXPDTC = as.Date(sample(
      c("2013-09-10", "2013-12-15", "2014-02-05", "2014-03-20"),
      dplyr::n(),
      replace = TRUE
    )),
    EXLOT = case_when(
      DAEXPDTC == as.Date("2013-09-10") ~ "XXX-YYY-ZZZ-004",
      DAEXPDTC == as.Date("2013-12-15") ~ "XXX-YYY-ZZZ-002",
      DAEXPDTC == as.Date("2014-02-05") ~ "XXX-YYY-ZZZ-003",
      DAEXPDTC == as.Date("2014-03-20") ~ "XXX-YYY-ZZZ-005",
      TRUE ~ "XXX-YYY-ZZZ-006" # Default value for all other dates
    ),
    ADOSE = case_when(
      is.na(EXDOSE) ~ NA_real_,
      EXDOSE == 54 ~ 7.5,
      EXDOSE == 81 ~ 15,
      .default = EXDOSE
    ),
    TRT01P = droplevels(dplyr::case_when(
      TRT01P == "Screen Failure" ~ NA,
      .default = TRT01P
    )),
    TRT01PN = dplyr::case_when(
      TRT01P == "Xanomeline High Dose" ~ 1,
      TRT01P == "Xanomeline Low Dose" ~ 2,
      TRT01P == "Placebo" ~ 3
    ),
    TRT01P = forcats::fct_reorder(TRT01P, TRT01PN, .na_rm = TRUE),
    TRT01A = droplevels(dplyr::case_when(
      TRT01A == "Screen Failure" ~ NA,
      .default = TRT01A
    )),
    AVISITN = case_when(
      toupper(VISIT) == "BASELINE" ~ 1,
      toupper(VISIT) == "WEEK 2" ~ 2,
      toupper(VISIT) == "WEEK 4" ~ 3,
      toupper(VISIT) == "WEEK 6" ~ 4,
      toupper(VISIT) == "WEEK 8" ~ 5,
      toupper(VISIT) == "WEEK 12" ~ 6,
      toupper(VISIT) == "WEEK 16" ~ 7,
      toupper(VISIT) == "WEEK 20" ~ 8,
      toupper(VISIT) == "WEEK 24" ~ 9,
      toupper(VISIT) == "WEEK 26" ~ 10,
    ),
    AVISIT = fct_reorder(
      as.factor(case_when(
        toupper(VISIT) == "BASELINE" ~ "Screening",
        toupper(VISIT) == "WEEK 2" ~ "Cycle 02",
        toupper(VISIT) == "WEEK 4" ~ "Cycle 03",
        toupper(VISIT) == "WEEK 6" ~ "Cycle 04",
        toupper(VISIT) == "WEEK 8" ~ "Cycle 05",
        toupper(VISIT) == "WEEK 12" ~ "Cycle 06",
        toupper(VISIT) == "WEEK 16" ~ "Cycle 07",
        toupper(VISIT) == "WEEK 20" ~ "Cycle 08",
        toupper(VISIT) == "WEEK 24" ~ "Cycle 09",
        toupper(VISIT) == "WEEK 26" ~ "End Of Treatment",
        TRUE ~ as.character(VISIT) # Other values remain unchanged
      )),
      AVISITN,
      .na_rm = FALSE
    ),
    TRT01AN = dplyr::case_when(
      TRT01A == "Xanomeline High Dose" ~ 1,
      TRT01A == "Xanomeline Low Dose" ~ 2,
      TRT01A == "Placebo" ~ 3
    ),
    TRT01A = forcats::fct_reorder(TRT01A, TRT01AN, .na_rm = TRUE),
    AGEGR1 = factor(
      dplyr::case_when(
        AGE >= 18 & AGE < 65 ~ ">=18 to <65",
        AGE >= 65 & AGE < 75 ~ ">=65 to <75",
        AGE >= 75 ~ ">=75"
      ),
      levels = c(
        ">=18 to <65",
        ">=65 to <75",
        ">=75"
      )
    ),
    AOCCUR = factor(
      dplyr::case_when(
        !is.na(EXDOSE) ~ "Y",
        .default = NA_character_
      ),
      levels = c("Y", "N")
    ),
    SEX = factor(
      dplyr::case_when(
        SEX == "F" ~ "Female",
        SEX == "M" ~ "Male"
      ),
      levels = c(
        "Female",
        "Male",
        "Intersex",
        "Unknown"
      )
    ),
    COUNTRY = as.factor("United States of America"),
    RACEGR1 = as.factor(RACEGR1),
    ETHNIC = factor(
      dplyr::case_when(
        ETHNIC == "HISPANIC OR LATINO" ~ "Hispanic or Latino",
        ETHNIC == "NOT HISPANIC OR LATINO" ~ "Not Hispanic or Latino",
        ETHNIC == "NOT REPORTED" ~ "Not reported",
        ETHNIC == "UNKNOWN" ~ "Unknown"
      ),
      levels = c(
        "Hispanic or Latino",
        "Not Hispanic or Latino",
        "Not reported",
        "Unknown"
      )
    ),
    ACAT1 = factor(
      case_when(
        is.na(EXADJ) ~ "Dose not adjusted",
        EXADJ == "ADVERSE EVENT" ~ "Dose adjusted",
        EXADJ == "MEDICAL ERROR" ~ "Dose adjusted",
        TRUE ~ "Other" # default case if needed
      ),
      levels = c(
        "Dose not adjusted",
        "Dose adjusted",
        "Other"
      )
    ),
    AREASOC = factor(
      case_when(
        is.na(EXADJ) ~ NA_character_,
        EXADJ == "ADVERSE EVENT" ~ "Adverse Event",
        EXADJ == "MEDICATION ERROR" ~ "Other",
        TRUE ~ NA_character_ # default case
      ),
      levels = c(
        "Adverse Event",
        "Other"
      )
    ),
    AREASOO = factor(
      case_when(
        AREASOC == "Other" ~ "Other reason",
        TRUE ~ NA_character_ # default case
      ),
      levels = c(
        "Other reason"
      )
    ),
    AACTDU = factor(
      sample(
        c(
          "INFUSION INTERRUPTED",
          NA_character_,
          "INFUSION RATE INCREASED",
          "INFUSION CONTINUED AT SAME RATE",
          "INFUSION ABORTED",
          "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      levels = c(
        "INFUSION INTERRUPTED",
        "INFUSION RATE INCREASED",
        "INFUSION CONTINUED AT SAME RATE",
        "INFUSION ABORTED",
        "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE"
      )
    ),
    AACTDU1 = factor(
      case_when(
        AACTDU == "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE" ~
          "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE",
        TRUE ~ NA_character_ # default case
      ),
      levels = c(
        "FULL DOSE ADMINISTERED WITHOUT INTERRUPTION OR RATE CHANGE"
      )
    ),
    AACTDU2 = case_when(
      AACTDU == "INFUSION ABORTED" ~ "INFUSION ABORTED",
      TRUE ~ NA_character_ # default case
    ),
    AACTDU3 = case_when(
      AACTDU == "INFUSION INTERRUPTED" ~ "INFUSION INTERRUPTED",
      TRUE ~ NA_character_ # default case
    ),
    AACTDU4 = case_when(
      AACTDU == "INFUSION RATE DECREASED" ~ "INFUSION RATE DECREASED",
      TRUE ~ NA_character_ # default case
    ),
    AACTDU5 = case_when(
      AACTDU == "INFUSION RATE INCREASED" ~ "INFUSION RATE INCREASED",
      TRUE ~ NA_character_ # default case
    ),
    ACAT2 = factor(
      sample(
        c("Dose not administered", "Dose administered"),
        dplyr::n(),
        replace = TRUE
      ),
      levels = c(
        "Dose not administered",
        "Dose administered"
      )
    ),
    AACTPR = factor(
      sample(
        c(
          "INFUSION RATE DECREASED COMPARED TO PRIOR INFUSION",
          NA_character_,
          "INFUSION SKIPPED (AND NOT MADE UP)",
          "DOSE RE-ESCALATED",
          "INFUSION DELAYED WITHIN THE CYCLE",
          "DOSE REDUCED COMPARED TO PRIOR INFUSION",
          "SAME DOSE AS PRIOR INFUSION",
          "STUDY DRUG PERMANENTLY DISCONTINUED",
          "DOSE REDUCED",
          "DOSE REDUCED COMPARED TO PRIOR INJECTION",
          "REDUCED [TARGET/TREATMENT]",
          "FULL DOSE ADMINISTERED WITH INTERRUPTION"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      levels = c(
        "INFUSION RATE DECREASED COMPARED TO PRIOR INFUSION",
        "INFUSION SKIPPED (AND NOT MADE UP)",
        "DOSE RE-ESCALATED",
        "INFUSION DELAYED WITHIN THE CYCLE",
        "DOSE REDUCED COMPARED TO PRIOR INFUSION",
        "SAME DOSE AS PRIOR INFUSION",
        "STUDY DRUG PERMANENTLY DISCONTINUED",
        "DOSE REDUCED",
        "DOSE REDUCED COMPARED TO PRIOR INJECTION",
        "REDUCED [TARGET/TREATMENT]",
        "FULL DOSE ADMINISTERED WITH INTERRUPTION"
      )
    ),
    AADJP = dplyr::case_when(
      !is.na(AACTPR) & as.character(AACTPR) != "" ~
        sample(
          c(levels(AREASOC), "Physician Decision", "Subject Request", "Administrative Reason", "Medical Reason"),
          dplyr::n(),
          replace = TRUE
        ),
      TRUE ~ NA_character_
    ),
    AADJPOTH = ifelse(
      !is.na(AADJP) & AADJP == "Other",
      sample(
        c(
          "Personal reasons",
          "Physician decision",
          "Subject request",
          "Administrative reasons",
          "Other medical reason"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      NA_character_
    ),
    AADJ = dplyr::case_when(
      !is.na(AACTDU) & as.character(AACTDU) != "" ~
        sample(
          c(levels(AREASOC), "Physician Decision", "Subject Request", "Administrative Reason", "Medical Reason"),
          dplyr::n(),
          replace = TRUE
        ),
      TRUE ~ NA_character_
    ),
    AADJOTH = ifelse(
      !is.na(AADJ) & AADJ == "Other",
      sample(
        c(
          "Cystoscopy with lesion. RTU pending",
          "Multiple tumors identified by cystoscopy on week 12.",
          "During week 12, the cystoscopy was performed and suspicion of progression was observed.",
          "Bladder mapping performed on 20 Nov 2023",
          "Progression",
          "Drug was never administered because the subject dropped out prior to first dose",
          "Personal reasons",
          "TAR200 was half blocked by the exit port of the UPC",
          "PQC",
          "The device could not push TAR-200 out.",
          "Sheet could not into the bladder",
          "Pt had his 24w TURB on 22th of December so TAR was removed at that time",
          "The dose delayed due to TURB performed 10-Nov-2023",
          "TAR-200 did not take the pretzel shape inside the bladder.",
          "Patient's vacation",
          "TAR200 was blocked in the middle of the UPC hole",
          "PQC was reported.",
          "PQC was identified",
          "Urinary placement catheter could not into the bladder.",
          "Insertion failed",
          "Suspicion of disease progression",
          "Due to personal schedule",
          "Urine culture not disponible"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      NA_character_
    ),
    ATPTN = dplyr::case_when(
      toupper(VISIT) == "BASELINE" ~ 0,
      toupper(VISIT) == "WEEK 2" ~ 1,
      toupper(VISIT) == "WEEK 4" ~ 2,
      toupper(VISIT) == "WEEK 6" ~ 3,
      toupper(VISIT) == "WEEK 8" ~ 4,
      toupper(VISIT) == "WEEK 12" ~ 5,
      toupper(VISIT) == "WEEK 16" ~ 6,
      toupper(VISIT) == "WEEK 20" ~ 7,
      toupper(VISIT) == "WEEK 24" ~ 8,
      toupper(VISIT) == "WEEK 26" ~ 9,
    ),
    ATPT = dplyr::case_when(
      ATPTN == 0 ~ "Morning dose",
      ATPTN == 1 ~ "Evening dose 1",
      ATPTN == 2 ~ "Evening dose 2",
      ATPTN == 3 ~ "Evening dose 3",
      ATPTN == 4 ~ "Evening dose 4",
      ATPTN == 5 ~ "Evening dose 5",
      ATPTN == 6 ~ "Evening dose 6",
      ATPTN == 7 ~ "Evening dose 7",
      ATPTN == 8 ~ "Evening dose 8",
      ATPTN == 9 ~ "Evening dose 9"
    ),
    ASCHDOSE = ifelse(!is.na(ECOCCUR), EXDOSE, NA_real_),
    ASCHDOSU = ifelse(!is.na(ASCHDOSE), EXDOSU, NA_character_),
    ADOSFRM = EXDOSFRM,
    ADOSU = EXDOSU,
    ADOSFRQ = EXDOSFRQ,
    ADOSFRQP = EXDOSFRQ,
    ASTDY = ASTDY,
    ECMOOD = "Scheduled",
    AROUTE = factor(
      sample(
        c(
          unique(na.omit(EXROUTE)),
          "Oral",
          "Injection",
          "Infusion"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      levels = c(
        unique(na.omit(EXROUTE)),
        "Oral",
        "Injection",
        "Infusion"
      )
    ),
    ATVINF = ADOSE,
    ATVINFU = ADOSU,
    AINFRAT = ADOSE,
    AINFRAU = ADOSU,
    EXSTDT = as.Date(EXSTDTC),
    # Add random hour times to ASTDTM
    ASTDTM = {
      # Extract the date part from ASTDTM
      date_parts <- sub(" UTC", "", ASTDTM)
      # Generate random hours (0-23) for each row
      random_hours <- sample(0:23, dplyr::n(), replace = TRUE)
      # Format the random hours as time strings (HH:00:00)
      random_times <- sprintf("%02d:00:00", random_hours)
      # Combine date with random time and UTC timezone
      paste(date_parts, random_times, "UTC")
    }
  )

  gen <- gen |>
    dplyr::mutate(
      ADOSU = ifelse(!is.na(ADOSE), EXDOSU, NA_character_)
    )

  # Derive APERIOD and APERIODC using overall (all subjects) min/max date midpoint
  .overall_min <- min(as.Date(sub(" .*", "", gen$ASTDTM)), na.rm = TRUE)
  .overall_max <- max(as.Date(sub(" .*", "", gen$ASTDTM)), na.rm = TRUE)
  .midpoint <- .overall_min + floor(as.numeric(.overall_max - .overall_min) / 2)

  gen <- gen |>
    dplyr::mutate(
      APERIOD = dplyr::if_else(
        as.Date(sub(" .*", "", ASTDTM)) <= .midpoint,
        1L,
        2L
      ),
      APERIODC = dplyr::if_else(APERIOD == 1L, "Period 1", "Period 2")
    )

  # TODO: Delete ------------
  # Start -------------------
  # Temporary Changes for ADISHUM Listings

  gen <- gen |>
    mutate(
      AVISITN = if_else(
        AVISITN == 9,
        3,
        AVISITN
      )
    )

  gen <- gen |>
    mutate(
      AVISIT = dplyr::case_when(
        AVISIT == "BASELINE" ~ "Day 1",
        AVISIT == "WEEK 2" ~ "Day 2",
        AVISIT == "WEEK 24" ~ "Day 3",
      )
    )

  gen <- gen |>
    mutate(
      AOCCUR = if_else(
        !is.na(EXSTDT),
        "Y",
        "N"
      )
    )

  # End -------------------

  # Derive ABODSYSy and ADECODy
  #  - ABODSYSy: AE System Organ Class (AE.AEBODSYS)
  #  - ADECODy:  AE Preferred Term    (AE.AEDECOD)

  # Actions indicating change due to AE (RELREC proxy)
  ae_actions <- c(
    "DRUG WITHDRAWN",
    "DRUG INTERRUPTED",
    "DOSE REDUCED",
    "DOSE RATE REDUCED"
  )

  ex_key <- gen |>
    dplyr::mutate(
      .row_id = dplyr::row_number(),
      EXASTDT = as.Date(sub(" .*", "", ASTDTM))
    ) |>
    dplyr::select(.row_id, USUBJID, EXASTDT)

  ae_sub <- adae |>
    dplyr::mutate(
      AEASTDT = as.Date(sub(" .*", "", ASTDTM))
    ) |>
    dplyr::filter(AEACN %in% ae_actions) |>
    dplyr::select(USUBJID, AEASTDT, AEBODSYS, AEDECOD, AESEQ)

  ex_ae <- dplyr::left_join(
    ex_key,
    ae_sub,
    by = c("USUBJID" = "USUBJID", "EXASTDT" = "AEASTDT")
  ) |>
    dplyr::arrange(.row_id, AESEQ)

  soc_summ <- ex_ae |>
    dplyr::group_by(.row_id) |>
    dplyr::summarise(
      ABODSYS1 = dplyr::first(unique(na.omit(AEBODSYS))),
      ABODSYS2 = dplyr::nth(unique(na.omit(AEBODSYS)), 2),
      ADECOD1 = dplyr::first(unique(na.omit(AEDECOD))),
      ADECOD2 = dplyr::nth(unique(na.omit(AEDECOD)), 2),
      .groups = "drop"
    )

  gen <- gen |>
    dplyr::mutate(.row_id = dplyr::row_number()) |>
    dplyr::left_join(soc_summ, by = ".row_id") |>
    dplyr::select(-.row_id)

  gen <- gen |>
    dplyr::mutate(
      AREASOC = dplyr::case_when(
        ABODSYS1 %in%
          c(
            "General disorders and administration site conditions",
            "Gastrointestinal disorders"
          ) ~ "Adverse Event",
        ABODSYS1 %in%
          c(
            "General disorders and administration site conditions",
            "Gastrointestinal disorders"
          ) ~ "Adverse Event",
        .default = AREASOC
      ),
      ADECOD1 = dplyr::case_when(
        toupper(AADJ) == "ADVERSE EVENT" | toupper(AADJP) == "ADVERSE EVENT" ~
          sample(
            c("VOMITING", "NAUSEA", "DIARRHOEA", "FATIGUE", "HEADACHE", "RASH"),
            dplyr::n(),
            replace = TRUE
          ),
        .default = NA_character_
      ),
      ADECOD2 = dplyr::case_when(
        toupper(AADJ) == "ADVERSE EVENT" | toupper(AADJP) == "ADVERSE EVENT" ~
          ifelse(
            sample(c(TRUE, FALSE), dplyr::n(), replace = TRUE, prob = c(0.5, 0.5)),
            sample(
              c("ANAEMIA", "NEUTROPENIA", "THROMBOCYTOPENIA", "PYREXIA", "OEDEMA PERIPHERAL"),
              dplyr::n(),
              replace = TRUE
            ),
            NA_character_
          ),
        .default = NA_character_
      ),
      ABODSYS1 = dplyr::if_else(
        !is.na(ADECOD1),
        sample(
          c(
            "Gastrointestinal disorders",
            "General disorders and administration site conditions",
            "Nervous system disorders",
            "Skin and subcutaneous tissue disorders"
          ),
          dplyr::n(),
          replace = TRUE
        ),
        NA_character_
      ),
      ABODSYS2 = dplyr::if_else(
        !is.na(ADECOD2),
        sample(
          c(
            "Investigations",
            "Musculoskeletal and connective tissue disorders",
            "Respiratory, thoracic and mediastinal disorders",
            "Vascular disorders"
          ),
          dplyr::n(),
          replace = TRUE
        ),
        NA_character_
      ),
      ANL01FL = dplyr::case_when(
        AROUTE == "Oral" & toupper(AACTPR) == "DOSE REDUCED" & ECMOOD == "Scheduled" ~ "Y",
        AROUTE == "Injection" &
          toupper(AACTPR) %in% c("DOSE REDUCED COMPARED TO PRIOR INJECTION", "REDUCED [TARGET/TREATMENT]") &
          ECMOOD == "Scheduled" ~ "Y",
        AROUTE == "Infusion" &
          toupper(AACTPR) %in% c("DOSE REDUCED COMPARED TO PRIOR INFUSION", "REDUCED [TARGET/TREATMENT]") &
          ECMOOD == "Scheduled" ~ "Y",
        .default = "N"
      ),
      ANL02FL = dplyr::case_when(
        ANL01FL == "Y" & AADJ == "Adverse Event" & ECMOOD == "Scheduled" ~ "Y",
        .default = "N"
      ),
      ANL03FL = dplyr::case_when(
        AROUTE == "Injection" & toupper(AACTPR) == "FULL DOSE ADMINISTERED WITH INTERRUPTION" ~ "Y",
        AROUTE == "Infusion" &
          (toupper(AACTDU) == "INFUSION INTERRUPTED" | toupper(AACTDU3) == "INFUSION INTERRUPTED") ~ "Y",
        .default = "N"
      ),
      ANL04FL = dplyr::case_when(
        AROUTE == "Injection" & toupper(AACTPR) == "INJECTION SKIPPED (AND NOT MADE UP)" ~ "Y",
        AROUTE == "Infusion" & toupper(AACTPR) == "INFUSION SKIPPED (AND NOT MADE UP)" ~ "Y",
        .default = "N"
      ),
      ACDOSE = ifelse(!is.na(ECOCCUR), ECDOSE, NA_real_),
      ACDOSU = factor(
        ifelse(!is.na(ACDOSE), "mL", NA_character_),
        levels = c("mL")
      ),
      ECRSDOSD = factor(
        sample(
          c(
            "Adverse Event",
            "Low Neutrophil Count",
            "Thrombocytopenia",
            "Anemia",
            "Infection",
            "Non-Hematologic Toxicity",
            "Recovery From Toxicity",
            "Pending Lab Results",
            "Physician Decision",
            "Subject Request",
            "Missed Visit",
            "Other"
          ),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c(
          "Adverse Event",
          "Low Neutrophil Count",
          "Thrombocytopenia",
          "Anemia",
          "Infection",
          "Non-Hematologic Toxicity",
          "Recovery From Toxicity",
          "Pending Lab Results",
          "Physician Decision",
          "Subject Request",
          "Missed Visit",
          "Other"
        )
      ),
      ECRSDSDO = ifelse(
        as.character(ECRSDOSD) == "Other",
        sample(
          c(
            "Weather conditions",
            "Site closure",
            "Equipment malfunction",
            "Staff unavailability"
          ),
          dplyr::n(),
          replace = TRUE
        ),
        NA_character_
      ),
      ADOSDLY = factor(
        sample(
          c("Y", "N", NA_character_),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c(
          "Y",
          "N"
        )
      ),
      ARSDOSD = ifelse(ADOSDLY == "Y" & ECMOOD == "Scheduled", as.character(ECRSDOSD), NA_character_),
      ARSDSDO = ifelse(ARSDOSD == "Other" & ECMOOD == "Scheduled", as.character(ECRSDSDO), NA_character_),
      ECAVAMT = factor(
        sample(
          c(1, 1.5, 1.7, 2.0, 5.0, 3.0, 4.0),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c(1, 1.5, 1.7, 2.0, 5.0, 3.0, 4.0)
      ),
      AVAMT = ifelse(!is.na(ECOCCUR), ECAVAMT, NA_real_),
      ECAVAMTU = factor(
        sample(
          c("mL"),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c("mL")
      ),
      AVAMTU = ifelse(!is.na(ECOCCUR), as.character(ECAVAMTU), NA_character_),
      ATDPRP = sample(c(100, 150, 200, 250, 300), n(), replace = TRUE),
      ATDPRPU = EXDOSU,
      ADURC = case_when(
        !is.na(AENDT) & !is.na(ASTDT) ~
          as.character(as.numeric(AENDT - ASTDT) + 1),
        !is.na(AENDTM) & !is.na(ASTDTM) ~
          sprintf(
            "%02d:%02d",
            as.numeric(difftime(AENDTM, ASTDTM, units = "hours")),
            as.numeric(difftime(AENDTM, ASTDTM, units = "mins")) %% 60
          ),
        TRUE ~ NA_character_
      )
    )

  gen <- gen |>
    dplyr::group_by(USUBJID, ATRT) |>
    dplyr::mutate(
      INTN = sum(ANL03FL == "Y", na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      INTCAT = dplyr::case_when(
        INTN == 1 ~ "1",
        INTN == 2 ~ "2",
        INTN >= 3 ~ ">=3",
        .default = NA_character_
      )
    )

  gen <- gen |>
    dplyr::group_by(USUBJID, ATRT) |>
    dplyr::mutate(
      SKPN = sum(ANL04FL == "Y", na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      SKPCAT = dplyr::case_when(
        SKPN == 1 ~ "1",
        SKPN == 2 ~ "2",
        SKPN >= 3 ~ ">=3",
        .default = NA_character_
      )
    )

  gen <- gen |>
    dplyr::group_by(USUBJID, ATRT) |>
    dplyr::mutate(
      DLYN = sum(ADOSDLY == "Y", na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      DLYCAT = dplyr::case_when(
        DLYN == 1 ~ "1",
        DLYN == 2 ~ "2",
        DLYN >= 3 ~ ">=3",
        .default = NA_character_
      )
    )

  # Join ADATRES and NABSTAT from ADISHUM - each collapsed separately, POSITIVE takes priority
  adatres_map <- adishum |>
    dplyr::filter(!is.na(ADATRES)) |>
    dplyr::mutate(.priority = dplyr::if_else(ADATRES == "POSITIVE", 1L, 2L)) |>
    dplyr::arrange(USUBJID, .priority) |>
    dplyr::distinct(USUBJID, .keep_all = TRUE) |>
    dplyr::select(USUBJID, ADATRES)

  nabstat_map <- adishum |>
    dplyr::filter(!is.na(NABSTAT)) |>
    dplyr::mutate(.priority = dplyr::if_else(NABSTAT == "POSITIVE", 1L, 2L)) |>
    dplyr::arrange(USUBJID, .priority) |>
    dplyr::distinct(USUBJID, .keep_all = TRUE) |>
    dplyr::select(USUBJID, NABSTAT)

  gen <- gen |>
    dplyr::left_join(adatres_map, by = "USUBJID") |>
    dplyr::left_join(nabstat_map, by = "USUBJID")

  # SECTION ---- Need specific value for EXSTDT for a Listing ------------
  # Get USUBJIDs that pass the ADISHUM filter and also exist in ADAE
  qualified_usubjids <- adishum |>
    dplyr::filter(
      IMFL == "Y",
      PARAMCD == "ADATRE",
      AVALC == "Y"
    ) |>
    dplyr::distinct(USUBJID) |>
    dplyr::filter(USUBJID %in% adae$USUBJID) |>
    dplyr::pull(USUBJID)

  # For each qualified subject, get one DOSEDT from ADAE and convert to Date
  # DOSEDT in ADAE is stored as integer (SAS date, origin 1960-01-01)
  # Result: one row per USUBJID with their DOSEDT_date
  dosedt_per_subj <- adae |>
    dplyr::filter(
      USUBJID %in% qualified_usubjids,
      !is.na(DOSEDT)
    ) |>
    dplyr::group_by(USUBJID) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      USUBJID,
      DOSEDT_date = as.Date(DOSEDT, origin = "1960-01-01")
    ) |>
    head(30)

  # Join DOSEDT_date onto ADEX by USUBJID, then set EXSTDT = DOSEDT_date
  # for the FIRST matching row per subject only (1 subject = 1 row updated)
  gen <- gen |>
    dplyr::left_join(dosedt_per_subj, by = "USUBJID") |>
    dplyr::group_by(USUBJID) |>
    dplyr::mutate(
      EXSTDT = dplyr::if_else(
        !is.na(DOSEDT_date) & dplyr::row_number() == 1,
        DOSEDT_date,
        EXSTDT
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-DOSEDT_date)

  # Additional labels for all relevant variables
  additional_labels <- list(
    EXTRT = "Planned Treatment",
    EXDOSE = "Adjusted Dose",
    ACAT1 = "Analysis Category 1",
    AADJOTH = "Other Anal Reason for Dose Adjustment",
    ACAT2 = "Analysis Category 2",
    AACTPR = "Action Taken Prior to Infusion Start",
    ANL01FL = "Analysis Flag 01-Dose reduction",
    ANL02FL = "Analysis Flag 02-Dose reductn due to AE",
    ANL03FL = "Analysis Flag 03-Dose interruption",
    ANL04FL = "Analysis Flag 04-Dose skip not made up",
    INTN = "Number of Dose Interruptions",
    INTCAT = "Dose Interruption Category",
    SKPN = "Number of Dose Skipped",
    SKPCAT = "Dose Skipped Category",
    AREASOC = "Analysis Reason for Occur Value",
    AADJ = "Analysis Reason for Dose Adjustment",
    AADJP = "Analysis Reason for Prior Dose Adjust",
    AREASOO = "Other Analysis Reason for Occur Value",
    AACTDU = "Analysis Action Taken During Study Trt",
    AACTDU1 = "Act Takn Dur Infus-Full Dose Admined",
    AACTDU2 = "Act Takn Dur Infus-Infusion Aborted",
    AACTDU3 = "Act Takn Dur Infus-Infusion Interrupted",
    AACTDU4 = "Act Takn Dur Infus-Infusion Rate Decrsed",
    AACTDU5 = "Act Takn Dur Infus-Infusion Rate Incrsed",
    TRT01P = "Planned Treatment for Period 01",
    TRT01A = "Actual Treatment for Period 01",
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    ARMCD = "Treatment Category Code",
    ARM = "Treatment Group",
    ACTARMCD = "Actual Arm Code",
    ACTARM = "Actual Treatment Group",
    DAEXPDTC = "Date of Exposure",
    EXLOT = "Lot Number",
    ADOSE = "Analysis Dose",
    AVISITN = "Visit Number",
    AVISIT = "Visit Label",
    AGEGR1 = "Age Group",
    AOCCUR = "Analysis Occurrence",
    SEX = "Sex",
    COUNTRY = "Country",
    RACE = "Race",
    ETHNIC = "Ethnicity",
    ATRT = "Analysis name of Treatment",
    ASCHDOSE = "Analysis Scheduled Dose",
    ASCHDOSU = "Analysis Scheduled Dose Units",
    ADOSFRM = "Analysis Dose Form",
    ADOSU = "Analysis Dose Units",
    ADOSFRQ = "Analysis Dosing Frequency per Interval",
    AROUTE = "Analysis Route of Administration",
    AADJPOTH = "Other Anal Reason for Dose Adjust Prior",
    ATVINF = "Analysis Total Volume Infused",
    ATVINFU = "Analysis Total Volume Infused Units",
    AINFRAT = "Analysis Infusion Rate",
    AINFRAU = "Analysis Infusion Rate Unit",
    ABODSYS1 = "AE SOC Driving Study Drug Action (1)",
    ABODSYS2 = "AE SOC Driving Study Drug Action (2)",
    ADECOD1 = "AE PT Driving Study Drug Action (1)",
    ADECOD2 = "AE PT Driving Study Drug Action (2)",
    ASTDY = "Analysis Study Day",
    ARSDOSD = "Analysis Reason for Dose Delayed",
    ARSDSDO = "Other Analysis Reason for Dose Delayed",
    ACDOSE = "Analysis Dose Collected",
    ACDOSU = "Analysis Dose Collected Units",
    ADOSDLY = "Analysis Dose Delayed",
    DLYN = "Number of Dose Delay",
    DLYCAT = "Dose Delay Category",
    ATPT = "Analysis Timepoint",
    ATPTN = "Analysis Timepoint (N)",
    ATDPRP = "Analysis Total Dose Prepared",
    ATDPRPU = "Analysis Total Dose Prepared Units",
    ECOCCUR = "Occurrence of Exposure",
    ADURC = "Analysis Duration (C)",
    ADOSFRQP = "Analysis Scheduled Dosing Frequency",
    AVAMT = "Analysis Injection Vol Prescribed",
    AVAMTU = "Analysis Injection Vol Prescribed Units",
    APERIOD = "Analysis Period",
    APERIODC = "Analysis Period (C)",
    ADATRES = "Treatment-emergent ADA Subject Status",
    NABSTAT = "NAB Status",
    EXSTDT = "Exposure Start Date"
  )

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
adex <- gen_adex()
