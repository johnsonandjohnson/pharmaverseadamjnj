# Load necessary libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(pharmaverseadam)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADDILI dataset
gen_addili <- function(seed = 123) {
  set.seed(seed)

  # Source ADSL and ADLB
  source(file.path("data-raw", "adlb.R"))
  source(file.path("data-raw", "adsl.R"))

  raw <- adlb

  # Define Variables to Keep
  adlb_vars <- c(
    "STUDYID", "USUBJID", "PARAMCD", "PARAM", "ADT", "ATM", "ADTM", "ADY",
    "AVISIT", "AVISITN", "ATPT", "ATPTN", "ATPTREF", "TRTP", "TRTPN",
    "TRTA", "TRTAN", "AVAL", "BASE", "ANRIND", "BNRIND", "ANRLO", "ANRHI",
    "ABLFL", "APOBLFL", "ONTRTFL", "LBSEQ", "VISIT", "VISITNUM", "LBNAM", "DTYPE"
  )

  adsl_vars <- c(
    "STUDYID", "USUBJID", "TRT01P", "TRT01PN", "TRT01A", "TRT01AN",
    "TRTSEQP", "TRTSEQPN", "TRTSEQA", "TRTSEQAN", "ARFSTDT", "ARFSTDTM",
    "TRTSDT", "TRTSDTM", "TRTEDT", "TRTEDTM", "AGE", "AGEU", "AGEGR1",
    "AGEGR1N", "SEX", "RACE", "COUNTRY", "SAFFL", "FASFL", "RANDFL",
    "DTHFL", "SITEID", "SUBJID"
  )

  # Modify subjects who have ALT and BILI records to have 2nd and 3rd quadrant data
  hys_subjs <- c("01-701-1028", "01-701-1034", "01-701-1047", "01-701-1115")
  chol_subj <- "01-701-1118"


  adlb <- raw |>
    group_by(USUBJID) |>
    mutate(
      .hys_alt_ast_mult = if_else(USUBJID %in% hys_subjs, rnorm(1, mean = 4.1, sd = 0.1), NA_real_),
      .hys_bili_mult = if_else(USUBJID %in% hys_subjs, rnorm(1, mean = 3.1, sd = 0.1), NA_real_),
      .chol_alp_mult = if_else(USUBJID == chol_subj, rnorm(1, mean = 2.5, sd = 0.08), NA_real_),
      .chol_bili_mult = if_else(USUBJID == chol_subj, rnorm(1, mean = 3.1, sd = 0.1), NA_real_),
      .chol_alt_ast_mult = if_else(USUBJID == chol_subj, rnorm(1, mean = 1.5, sd = 0.05), NA_real_),
      AVAL = case_when(
        # Hy's Law: ALT & AST > 3x, BILI > 2x
        USUBJID %in% hys_subjs & PARAMCD %in% c("ALT", "AST") ~ ANRHI * .hys_alt_ast_mult,
        USUBJID %in% hys_subjs & PARAMCD == "BILI" ~ ANRHI * .hys_bili_mult,
        # Cholestasis: ALP > 2x, BILI > 2x
        USUBJID == chol_subj & PARAMCD == "ALP" ~ ANRHI * .chol_alp_mult,
        USUBJID == chol_subj & PARAMCD == "BILI" ~ ANRHI * .chol_bili_mult,
        USUBJID == chol_subj & PARAMCD %in% c("ALT", "AST") ~ ANRHI * .chol_alt_ast_mult,
        TRUE ~ AVAL
      ),
      R2ANRHI = AVAL / ANRHI
    ) |>
    ungroup() |>
    select(-.hys_alt_ast_mult, -.hys_bili_mult, -.chol_alp_mult, -.chol_bili_mult, -.chol_alt_ast_mult)

  # Subset ADLB and Keep Predecessor Variables
  addili_base <- adlb |>
    filter(
      PARAMCD %in% c("ALP", "ALT", "AST", "BILI"),
      !is.na(AVAL)
    ) |>
    select(any_of(adlb_vars))

  # join ADSL
  adsl_subset <- adsl |>
    select(any_of(adsl_vars))

  addili_merged <- addili_base |>
    left_join(
      adsl_subset,
      by = c("STUDYID", "USUBJID")
    )

  addili_prep <- addili_merged |>
    mutate(
      PARCAT1 = case_when(
        PARAMCD %in% c("ALT", "AST") ~ "ALT or AST",
        PARAMCD == "BILI" ~ "TBILI",
        PARAMCD == "ALP" ~ "ALP",
        TRUE ~ as.character(PARAMCD)
      ),
      R2ANRHI = AVAL / ANRHI,
      CRIT1 = case_when(
        PARCAT1 == "ALT or AST" ~ ">= 3x ULN",
        PARCAT1 %in% c("TBILI", "ALP") ~ ">= 2x ULN",
        TRUE ~ NA_character_
      ),
      CRIT1FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 3 ~ "Y",
        PARCAT1 %in% c("TBILI", "ALP") & R2ANRHI >= 2 ~ "Y",
        !is.na(R2ANRHI) & !is.na(CRIT1) ~ "N",
        TRUE ~ NA_character_
      ),
      CRIT2 = case_when(
        PARCAT1 == "ALT or AST" ~ ">= 5x ULN",
        PARCAT1 == "ALP" ~ ">= 1.5x ULN",
        PARCAT1 == "BILI" ~ ">= 3x ULN",
        TRUE ~ NA_character_
      ),
      CRIT2FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 5 ~ "Y",
        PARCAT1 == "ALT or AST" & !is.na(R2ANRHI) ~ "N",
        PARCAT1 == "ALP" & R2ANRHI >= 1.5 ~ "Y",
        PARCAT1 == "BILI" & R2ANRHI >= 3 ~ "Y",
        !is.na(R2ANRHI) & !is.na(CRIT2) ~ "N",
        TRUE ~ NA_character_
      ),
      CRIT3 = case_when(
        PARCAT1 == "ALT or AST" ~ ">= 8x ULN",
        PARCAT1 == "TBILI" ~ ">= 5x ULN",
        PARCAT1 == "ALP" ~ ">= 3x ULN",
        TRUE ~ NA_character_
      ),
      CRIT3FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 8 ~ "Y",
        PARCAT1 == "ALT or AST" & !is.na(R2ANRHI) ~ "N",
        PARCAT1 == "TBILI" & R2ANRHI >= 5 ~ "Y",
        PARCAT1 == "ALP" & R2ANRHI >= 3 ~ "Y",
        !is.na(R2ANRHI) & !is.na(CRIT3) ~ "N",
        TRUE ~ NA_character_
      ),
      CRIT4 = if_else(PARCAT1 == "ALT or AST", ">= 3x ULN and <5x ULN", NA_character_),
      CRIT4FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 3 & R2ANRHI < 5 ~ "Y",
        PARCAT1 == "ALT or AST" & !is.na(R2ANRHI) ~ "N",
        TRUE ~ NA_character_
      ),
      CRIT5 = if_else(PARCAT1 == "ALT or AST", ">= 5x ULN and <8x ULN", NA_character_),
      CRIT5FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 5 & R2ANRHI < 8 ~ "Y",
        PARCAT1 == "ALT or AST" & !is.na(R2ANRHI) ~ "N",
        TRUE ~ NA_character_
      ),
      CRIT6 = if_else(PARCAT1 == "ALT or AST", ">= 10x ULN", NA_character_),
      CRIT6FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 10 ~ "Y",
        PARCAT1 == "ALT or AST" & !is.na(R2ANRHI) ~ "N",
        TRUE ~ NA_character_
      ),
      CRIT7 = if_else(PARCAT1 == "ALT or AST", ">= 20x ULN", NA_character_),
      CRIT7FL = case_when(
        PARCAT1 == "ALT or AST" & R2ANRHI >= 20 ~ "Y",
        PARCAT1 == "ALT or AST" & !is.na(R2ANRHI) ~ "N",
        TRUE ~ NA_character_
      )
    )

  # Derive ANL04FL, ANL06FL
  anl04 <- addili_prep |>
    filter(PARCAT1 == "ALT or AST", ONTRTFL == "Y", CRIT1FL == "Y") |>
    pull(USUBJID) |>
    unique()

  anl06 <- addili_prep |>
    filter(PARCAT1 == "ALP", ONTRTFL == "Y", CRIT1FL == "Y") |>
    pull(USUBJID) |>
    unique()

  addili_anl0406 <- addili_prep |>
    mutate(
      ANL04FL = if_else(USUBJID %in% anl04, "Y", NA_character_),
      ANL06FL = if_else(USUBJID %in% anl06, "Y", NA_character_)
    )
  # Derive ANL05FL, ANL07FL
  anl05 <- addili_anl0406 |>
    filter(PARCAT1 == "ALT or AST", ONTRTFL == "Y", CRIT1FL == "Y") |>
    select(USUBJID, TRG_ADT = ADT)

  anl07 <- addili_anl0406 |>
    filter(PARCAT1 == "ALP", ONTRTFL == "Y", CRIT1FL == "Y") |>
    select(USUBJID, TRG_ADT = ADT)

  # Apply HDILI temporal window (30 days here)
  addili_anl04050607 <- derive_var_joined_exist_flag(
    addili_anl0406,
    dataset_add = anl05,
    by_vars = exprs(USUBJID),
    new_var = ANL05_CHECK,
    join_vars = exprs(TRG_ADT),
    join_type = "all",
    filter_join = ADT >= TRG_ADT & ADT <= TRG_ADT + 30
  ) |>
    mutate(
      ANL05FL = if_else(PARAMCD == "BILI" & ANL04FL == "Y" & ANL05_CHECK == "Y", "Y", NA_character_)
    )

  # Apply CDILI temporal window (30 days here)
  addili_anl04050607 <- derive_var_joined_exist_flag(
    addili_anl04050607,
    dataset_add = anl07,
    by_vars = exprs(USUBJID),
    new_var = ANL07_CHECK,
    join_vars = exprs(TRG_ADT),
    join_type = "all",
    filter_join = ADT >= TRG_ADT & ADT <= TRG_ADT + 30
  ) |>
    mutate(
      ANL07FL = if_else(PARAMCD == "BILI" & ANL06FL == "Y" & ANL07_CHECK == "Y", "Y", NA_character_)
    ) |>
    select(-ANL05_CHECK, -ANL07_CHECK) # remove temporary variables


  # Derive ANL02FL - Max R2ANRHI for HDILI
  addili_anl0204050607 <- restrict_derivation(
    addili_anl04050607,
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARCAT1),
      order = exprs(desc(R2ANRHI), ADT, LBSEQ, AVISITN),
      new_var = ANL02FL,
      mode = "first",
      check_type = "none"
    ),
    filter = (PARAMCD == "BILI" & (ONTRTFL == "Y" | ANL05FL == "Y")) |
      (PARAMCD != "BILI" & ONTRTFL == "Y")
  )

  # Derive ANL03FL - Max R2ANRHI for CDILI
  addili_anl020304050607 <- restrict_derivation(
    addili_anl0204050607,
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARCAT1),
      order = exprs(desc(R2ANRHI), ADT, LBSEQ, AVISITN),
      new_var = ANL03FL,
      mode = "first",
      check_type = "none"
    ),
    filter = (PARAMCD == "BILI" & (ONTRTFL == "Y" | ANL07FL == "Y")) |
      (PARAMCD != "BILI" & ONTRTFL == "Y")
  )

  subj_summary <- addili_anl020304050607 |>
    group_by(STUDYID, USUBJID) |>
    summarise(
      across(any_of(c(
        "TRT01P", "TRT01PN", "TRT01A", "TRT01AN", "AGE", "AGEU", "AGEGR1", "AGEGR1N", "SEX",
        "RACE", "COUNTRY", "SAFFL", "FASFL", "RANDFL", "SITEID", "SUBJID"
      )), first),
      has_ontrt_alt_ast = any(PARCAT1 == "ALT or AST" & ONTRTFL == "Y", na.rm = TRUE),
      has_elig_bili_hdili = any(PARAMCD == "BILI" & (ONTRTFL == "Y" | ANL05FL == "Y"), na.rm = TRUE),
      has_ontrt_alp = any(PARCAT1 == "ALP" & ONTRTFL == "Y", na.rm = TRUE),
      has_elig_bili_cdili = any(PARAMCD == "BILI" & (ONTRTFL == "Y" | ANL07FL == "Y"), na.rm = TRUE),

      # Create subject-level analysis flags
      ANL04FL = if_else(any(ANL04FL == "Y", na.rm = TRUE), "Y", NA_character_),
      ANL06FL = if_else(any(ANL06FL == "Y", na.rm = TRUE), "Y", NA_character_),
      has_anl05fl = any(ANL05FL == "Y", na.rm = TRUE),
      has_anl07fl = any(ANL07FL == "Y", na.rm = TRUE),

      #  BILI CRIT1 rules
      bili_crit1_hdili = any(PARAMCD == "BILI" & CRIT1FL == "Y" & (ONTRTFL == "Y" | ANL05FL == "Y"), na.rm = TRUE),
      bili_crit1_cdili = any(PARAMCD == "BILI" & CRIT1FL == "Y" & (ONTRTFL == "Y" | ANL07FL == "Y"), na.rm = TRUE),
      bili_crit1_ontrt = any(PARAMCD == "BILI" & CRIT1FL == "Y" & ONTRTFL == "Y", na.rm = TRUE),
      .groups = "drop"
    )

  # ==============================================================================
  # Step 10: Derive HDILI Virtual Records + AVALCATy Mapping
  # ==============================================================================
  hdili_records <- subj_summary |>
    filter(has_ontrt_alt_ast & has_elig_bili_hdili) |>
    mutate(
      PARAMCD = "HDILI",
      PARAM = "Analysis Value - HDILI",
      PARCAT1 = NA_character_,
      ANL07FL = NA_character_,
      AVALC = factor(
        case_when(
          bili_crit1_hdili & !is.na(ANL04FL) ~ "ALT or AST >=3x ULN and TBILI >=2x ULN",
          bili_crit1_ontrt & is.na(ANL04FL) ~ "ALT and AST <3x ULN and TBILI >=2x ULN",
          !is.na(ANL04FL) ~ "ALT or AST >=3x ULN and TBILI <2x ULN",
          TRUE ~ "ALT and AST <3x ULN and TBILI <2x ULN"
        ),
        levels = c(
          "ALT or AST >=3x ULN and TBILI >=2x ULN",
          "ALT and AST <3x ULN and TBILI >=2x ULN",
          "ALT or AST >=3x ULN and TBILI <2x ULN",
          "ALT and AST <3x ULN and TBILI <2x ULN"
        )
      ),
      CRIT1 = "TBILI >=2x ULN within 30 days on or after ALT or AST >=3x ULN",
      CRIT1FL = if_else(AVALC == "ALT or AST >=3x ULN and TBILI >=2x ULN" & has_anl05fl, "Y", "N")
    )

  # map the eDISH quadrants for HDILI
  hdili_def <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N, ~AVALCAT2, ~AVALCA2N,
    AVALC == "ALT or AST >=3x ULN and TBILI >=2x ULN", "Potential Hy's Law", 1, "right upper quadrant", 1,
    AVALC == "ALT and AST <3x ULN and TBILI >=2x ULN", "Cholestasis", 2, "left upper quadrant", 2,
    AVALC == "ALT or AST >=3x ULN and TBILI <2x ULN", "Temple's Corollary", 3, "right lower quadrant", 3,
    AVALC == "ALT and AST <3x ULN and TBILI <2x ULN", "Normal or near normal", 4, "left lower quadrant", 4
  )
  hdili_records <- derive_vars_cat(hdili_records, definition = hdili_def)

  cdili_records <- subj_summary |>
    filter(has_ontrt_alp & has_elig_bili_cdili) |>
    mutate(
      PARAMCD = "CDILI",
      PARAM = "Analysis Value - CDILI",
      PARCAT1 = NA_character_,
      AVALC = case_when(
        bili_crit1_cdili & !is.na(ANL06FL) ~ "TBILI >=2x ULN and ALP >=2x ULN",
        bili_crit1_ontrt & is.na(ANL06FL) ~ "TBILI >=2x ULN and ALP <2x ULN",
        !is.na(ANL06FL) ~ "TBILI <2x ULN and ALP >=2x ULN",
        TRUE ~ "TBILI <2x ULN and ALP <2x ULN"
      ),
      CRIT1 = "TBILI >=2x ULN within 30 days on or after ALP >=2x ULN",
      CRIT1FL = if_else(AVALC == "TBILI >=2x ULN and ALP >=2x ULN" & has_anl07fl, "Y", "N")
    )

  cdili_def <- exprs(
    ~condition, ~AVALCAT2, ~AVALCA2N,
    AVALC == "TBILI >=2x ULN and ALP >=2x ULN", "right upper quadrant", 1,
    AVALC == "TBILI >=2x ULN and ALP <2x ULN", "left upper quadrant", 2,
    AVALC == "TBILI <2x ULN and ALP >=2x ULN", "right lower quadrant", 3,
    AVALC == "TBILI <2x ULN and ALP <2x ULN", "left lower quadrant", 4
  )
  cdili_records <- derive_vars_cat(cdili_records, definition = cdili_def)


  gen <- bind_rows(
    addili_anl020304050607,
    hdili_records |> select(
      any_of(names(addili_anl020304050607)),
      AVALC, AVALCAT1, AVALCA1N, AVALCAT2, AVALCA2N, CRIT1, CRIT1FL
    ),
    cdili_records |> select(
      any_of(names(addili_anl020304050607)),
      AVALC, AVALCAT2, AVALCA2N, CRIT1, CRIT1FL
    )
  ) |>
    # Derive ASEQ
    derive_var_obs_number(
      new_var = ASEQ,
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(PARAMCD, AVISITN, ADT)
    ) |>
    arrange(STUDYID, USUBJID, PARAMCD, AVISITN, ADT)


  additional_labels <- list(
    STUDYID = "Study Identifier",
    USUBJID = "Unique Subject Identifier",
    PARAMCD = "Parameter Code",
    PARAM = "Parameter",
    ADT = "Analysis Date",
    ADTM = "Analysis Date/Time",
    ADY = "Analysis Relative Day",
    AVISIT = "Analysis Visit",
    AVISITN = "Analysis Visit (N)",
    ATPT = "Analysis Timepoint",
    TRTP = "Planned Treatment",
    TRTA = "Actual Treatment",
    AVAL = "Analysis Value",
    BASE = "Baseline Value",
    ANRIND = "Analysis Reference Range Indicator",
    BNRIND = "Baseline Reference Range Indicator",
    ANRLO = "Analysis Normal Range Lower Limit",
    ANRHI = "Analysis Normal Range Upper Limit",
    ABLFL = "Baseline Record Flag",
    APOBLFL = "Post-Baseline Record Flag",
    ONTRTFL = "On Treatment Record Flag",
    LBSEQ = "Sequence Number",
    VISIT = "Visit Name",
    VISITNUM = "Visit Number",
    LBNAM = "Laboratory Name",
    DTYPE = "Derivation Type",
    TRT01P = "Planned Treatment for Period 01",
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01A = "Actual Treatment for Period 01",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    TRTSDT = "Date of First Exposure to Treatment",
    TRTSDTM = "Datetime of First Exposure to Treatment",
    TRTEDT = "Date of Last Exposure to Treatment",
    TRTEDTM = "Datetime of Last Exposure to Treatment",
    AGE = "Age",
    AGEU = "Age Units",
    AGEGR1 = "Pooled Age Group 1",
    AGEGR1N = "Pooled Age Group 1 (N)",
    SEX = "Sex",
    RACE = "Race",
    COUNTRY = "Country",
    SAFFL = "Safety Population Flag",
    FASFL = "Full Analysis Set Population Flag",
    RANDFL = "Randomized Population Flag",
    DTHFL = "Subject Death Flag",
    SITEID = "Study Site Identifier",
    SUBJID = "Subject Identifier for the Study",
    PARCAT1 = "Parameter Category 1",
    R2ANRHI = "Ratio to ANR Upper Limit",
    CRIT1 = "Analysis Criterion 1",
    CRIT1FL = "Criterion 1 Evaluation Result Flag",
    CRIT2 = "Analysis Criterion 2",
    CRIT2FL = "Criterion 2 Evaluation Result Flag",
    CRIT3 = "Analysis Criterion 3",
    CRIT3FL = "Criterion 3 Evaluation Result Flag",
    CRIT4 = "Analysis Criterion 4",
    CRIT4FL = "Criterion 4 Evaluation Result Flag",
    CRIT5 = "Analysis Criterion 5",
    CRIT5FL = "Criterion 5 Evaluation Result Flag",
    CRIT6 = "Analysis Criterion 6",
    CRIT6FL = "Criterion 6 Evaluation Result Flag",
    CRIT7 = "Analysis Criterion 7",
    CRIT7FL = "Criterion 7 Evaluation Result Flag",
    ANL04FL = "Analysis 04 Record Flag",
    ANL06FL = "Analysis 06 Record Flag",
    ANL05FL = "Analysis 05 Record Flag",
    ANL07FL = "Analysis 07 Record Flag",
    ANL02FL = "Analysis 02 Record Flag",
    ANL03FL = "Analysis 03 Record Flag",
    AVALC = "Analysis Value (C)",
    AVALCAT1 = "Analysis Value Category 1",
    AVALCA1N = "Analysis Value Category 1 (N)",
    AVALCAT2 = "Analysis Value Category 2",
    AVALCA2N = "Analysis Value Category 2 (N)",
    ASEQ = "Analysis Sequence Number"
  )

  # Standardize all Criterion flags to factors with Y first
  gen <- gen |>
    mutate(across(
      starts_with("CRIT") & ends_with("FL"),
      ~ factor(.x, levels = c("Y", "N"))
    ))

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

addili <- gen_addili()
