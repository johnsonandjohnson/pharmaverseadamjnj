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
      is.na(DOSEO) ~ NA_real_,
      DOSEO == 0 ~ 60,
      DOSEO < 60 ~ 60,
      DOSEO >= 60 & DOSEO < 180 ~ 120,
      DOSEO >= 180 & DOSEO < 240 ~ 180,
      DOSEO >= 240 & DOSEO < 300 ~ 240,
      DOSEO >= 300 & DOSEO < 480 ~ 300,
      DOSEO >= 480 & DOSEO < 720 ~ 480,
      DOSEO >= 720 & DOSEO < 1000 ~ 720,
      DOSEO >= 1000 & DOSEO < 5000 ~ 960,
      DOSEO >= 5000 & DOSEO < 10000 ~ 1800,
      DOSEO >= 10000 ~ 3600
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
    AOCCUR = as.factor(sample(c("N", "Y"), dplyr::n(), replace = TRUE)),
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
    AADJ = AREASOC,
    AADJPOTH = dplyr::case_when(
      AADJ == "Other" ~ "************",
      .default = "Reason prior to infusion"
    ),
    AADJP = AADJ,
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
    AADJOTH = ifelse(
      AREASOC == "Other",
      sample(
        c(
          "CYSTOSCOPI WITH LESION. RTU PENDING",
          "MULTIPLE TUMORS IDENTIFIED BY CYSTOSCOPY ON WEEK 12.",
          "DURING WEEK 12, THE CYSTOSCOPY WAS PERFORMED AND SUSPICION OF PROGRESSION WAS OBSERVED.",
          "BLADDER MAPPING PERFORMED ON 20 NOV 2023",
          "PROGRESSION",
          "DRUG WAS NEVER ADMINISTERED BECAUSE THE SUBJECT DROPPED OUT PRIOR TO FIRST DOSE",
          "PERSONAL REASONS",
          "TAR200 WAS HALF BLOCKED BY THE EXIT PORT OF THE UPC",
          "PQC",
          "THE DEVICE COULD NOT PUSH TAR-200 OUT.",
          "SHEET COULD NOT INTO THE BLADDER",
          "PT HAD HIS 24W TURB ON 22TH OF DECEMBER SO TAR WAS REMOVED AT THAT TIME",
          "THE DOSE DELAYED DUE TO TURB PERFORMED 10-NOV-2023",
          "TAR-200 DID NOT TAKE THE PRETZEL SHAPE INSIDE THE BLADDER.",
          "PATIENT'S VACATION",
          "TAR200 WAS BLOCKED IN THE MIDDLE OF THE UPC HOLE",
          "PQC WAS REPORTED.",
          "PQC WAS IDENTIFIED",
          "URINARY PLACEMENT CATHETER COULD NOT INTO THE BLADDER.",
          "INSERTION FAILED",
          "SUSPICION OF DISEASE PROGRESSION",
          "DUE TO PERSONAL SCHEDULE",
          "URINE CULTURE NOT DISPONIBLE"
        ),
        dplyr::n(),
        replace = TRUE
      ),
      "Reason during infusion"
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
      ATPTN == 0 ~ "Pre-dose",
      ATPTN == 1 ~ "Post-dose Hour 1",
      ATPTN == 2 ~ "Post-dose Hour 2",
      ATPTN == 3 ~ "Post-dose Hour 3",
      ATPTN == 4 ~ "Post-dose Hour 4",
      ATPTN == 5 ~ "Post-dose Hour 5",
      ATPTN == 6 ~ "Post-dose Hour 6",
      ATPTN == 7 ~ "Post-dose Hour 7",
      ATPTN == 8 ~ "Post-dose Hour 8",
      ATPTN == 9 ~ "Post-dose Hour 9"
    ),
    ASCHDOSE = EXDOSE,
    ASCHDOSU = EXDOSU,
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
        .default = AREASOC
      ),
      ABODSYS1 = dplyr::case_when(
        AADJ == "Adverse Event" ~ "Gastrointestinal disorders",
        .default = ABODSYS1
      ),
      ADECOD1 = dplyr::case_when(
        AADJ == "Adverse Event" ~ "VOMITING",
        .default = ADECOD1
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
      ECOCCUR = factor(
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
      ECDOSE = factor(
        sample(
          c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
      ),
      ACDOSE = ifelse(!is.na(ECOCCUR), ECDOSE, NA_real_),
      ACDOSU = factor(
        sample(
          c("mL"),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c("mL")
      ),
      ECRSDOSD = factor(
        sample(
          c(
            "ADVERSE EVENT",
            "LOW NEUTROPHIL COUNT",
            "THROMBOCYTOPENIA",
            "ANEMIA",
            "INFECTION",
            "NON-HEMATOLOGIC TOXICITY",
            "RECOVERY FROM TOXICITY",
            "PENDING LAB RESULTS",
            "PHYSICIAN DECISION",
            "SUBJECT REQUEST",
            "MISSED VISIT",
            "OTHER"
          ),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c(
          "ADVERSE EVENT",
          "LOW NEUTROPHIL COUNT",
          "THROMBOCYTOPENIA",
          "ANEMIA",
          "INFECTION",
          "NON-HEMATOLOGIC TOXICITY",
          "RECOVERY FROM TOXICITY",
          "PENDING LAB RESULTS",
          "PHYSICIAN DECISION",
          "SUBJECT REQUEST",
          "MISSED VISIT",
          "OTHER"
        )
      ),
      ECRSDSDO = ifelse(
        toupper(as.character(ECRSDOSD)) == "OTHER",
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
      ARSDOSD = ifelse(ECMOOD == "Scheduled", as.character(ECRSDOSD), NA_character_),
      ARSDSDO = ifelse(ECMOOD == "Scheduled", as.character(ECRSDSDO), NA_character_),
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
      ECAVAMT = factor(
        sample(
          c(1, 1.5, 1.7, 2.0, 5.0, 3.0, 4.0),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c(1, 1.5, 1.7, 2.0, 5.0, 3.0, 4.0)
      ),
      AVAMT = ifelse(is.na(ECOCCUR), ECAVAMT, NA_real_),
      ECAVAMTU = factor(
        sample(
          c("mL"),
          dplyr::n(),
          replace = TRUE
        ),
        levels = c("mL")
      ),
      AVAMTU = ifelse(is.na(ECOCCUR), as.character(ECAVAMTU), NA_character_),
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
    ACDOSE = "Analysis Dose Collected",
    ACDOSU = "Analysis Dose Collected Units",
    ADOSDLY = "Analysis Dose Delayed",
    DLYN = "Number of Dose Delay",
    DLYCAT = "Dose Delay Category",
    ATPT = "Analysis Timepoint",
    ATPTN = "Analysis Timepoint (N)",
    ATDPRP = "Analysis Total Dose Prepared",
    ATDPRPU = "Analysis Total Dose Prepared Units",
    ADURC = "Analysis Duration (C)",
    ADOSFRQP = "Analysis Scheduled Dosing Frequency",
    AVAMT = "Analysis Injection Vol Prescribed",
    AVAMTU = "Analysis Injection Vol Prescribed Units"
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
