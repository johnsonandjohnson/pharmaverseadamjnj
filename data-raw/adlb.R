# Generate ADLB dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(formatters)
library(forcats)
library(admiral)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADLB dataset
gen_adlb <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  # Get source data
  raw <- pharmaverseadam::adlb
  gen <- dplyr::select(raw, -AGEGR1)

  gen <- dplyr::mutate(
    gen,
    # We'll adjust AVAL and ANRHI for ALKPH to ensure ratio > 3 when calculated on-the-fly
    AVAL = ifelse(
      PARAMCD == "ALKPH" & USUBJID == unique(USUBJID)[1] & ADY > 0,
      30, # Use a high AVAL
      AVAL
    ),
    ANRHI = ifelse(
      PARAMCD == "ALKPH" & USUBJID == unique(USUBJID)[1] & ADY > 0,
      5, # This will ensure AVAL/ANRHI = 30/5 = 6 > 3
      ANRHI
    ),
    # Treatment and arm variables
    TRT01P = forcats::fct_recode(
      TRT01P,
      "Apalutamide" = "Xanomeline High Dose",
      "Apalutamide Subgroup" = "Xanomeline Low Dose"
    ),
    ARMCD = as.factor(dplyr::case_when(
      ARMCD == "Xan_Hi" ~ "Apa",
      ARMCD == "Xan_Lo" ~ "Apa_Sub",
      .default = ARMCD
    )),
    ARM = as.factor(dplyr::case_when(
      ARM == "Xanomeline High Dose" ~ "Apalutamide",
      ARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
      .default = ARM
    )),
    ACTARMCD = as.factor(dplyr::case_when(
      ACTARMCD == "Xan_Hi" ~ "Apa",
      ACTARMCD == "Xan_Lo" ~ "Apa_Sub",
      .default = ACTARMCD
    )),
    ACTARM = as.factor(dplyr::case_when(
      ACTARM == "Xanomeline High Dose" ~ "Apalutamide",
      ACTARM == "Xanomeline Low Dose" ~ "Apalutamide Subgroup",
      .default = ACTARM
    )),
    TRT01P = droplevels(dplyr::case_when(
      TRT01P == "Screen Failure" ~ NA,
      .default = TRT01P
    )),
    TRT01PN = dplyr::case_when(
      TRT01P == "Apalutamide" ~ 1,
      TRT01P == "Apalutamide Subgroup" ~ 2,
      TRT01P == "Placebo" ~ 3
    ),
    TRT01P = forcats::fct_reorder(TRT01P, TRT01PN, .na_rm = TRUE),
    TRT01A = forcats::fct_recode(
      TRT01A,
      "Apalutamide" = "Xanomeline High Dose",
      "Apalutamide Subgroup" = "Xanomeline Low Dose"
    ),
    TRT01A = droplevels(dplyr::case_when(
      TRT01A == "Screen Failure" ~ NA,
      .default = TRT01A
    )),
    TRT01AN = dplyr::case_when(
      TRT01A == "Apalutamide" ~ 1,
      TRT01A == "Apalutamide Subgroup" ~ 2,
      TRT01A == "Placebo" ~ 3
    ),

    # Analysis values
    AVAL = round(AVAL, 4),
    AVALC = as.character(AVAL),
    AVALU = LBSTRESU,
    ANL02FL = "Y",

    # Visit variables
    AVISITN = case_when(
      AVISIT == "Baseline" ~ 1,
      AVISIT == "Week 2" ~ 2,
      AVISIT == "Week 4" ~ 3,
      AVISIT == "Week 6" ~ 4,
      AVISIT == "Week 8" ~ 5,
      AVISIT == "Unscheduled 5.1" ~ 6,
      AVISIT == "Week 12" ~ 7,
      AVISIT == "Unscheduled 6.1" ~ 8,
      AVISIT == "Week 16" ~ 9,
      AVISIT == "Unscheduled 7.1" ~ 10,
      AVISIT == "Week 20" ~ 11,
      AVISIT == "Unscheduled 8.2" ~ 12,
      AVISIT == "Week 24" ~ 13,
      AVISIT == "Unscheduled 9.2" ~ 14,
      AVISIT == "Unscheduled 9.3" ~ 15,
      AVISIT == "Unscheduled 12.1" ~ 16,
      AVISIT == "Unscheduled 13.1" ~ 17,
      AVISIT == "Week 26" ~ 18
    ),
    AVISIT = fct_reorder(
      as.factor(case_when(
        AVISIT == "Baseline" ~ "Baseline",
        AVISIT == "Week 2" ~ "Cycle 02",
        AVISIT == "Week 4" ~ "Cycle 03",
        AVISIT == "Week 6" ~ "Cycle 04",
        AVISIT == "Week 8" ~ "Cycle 05",
        AVISIT == "Unscheduled 5.1" ~ "Cycle 06",
        AVISIT == "Week 12" ~ "Cycle 07",
        AVISIT == "Unscheduled 6.1" ~ "Cycle 08",
        AVISIT == "Week 16" ~ "Cycle 09",
        AVISIT == "Unscheduled 7.1" ~ "Cycle 10",
        AVISIT == "Week 20" ~ "Cycle 11",
        AVISIT == "Unscheduled 8.2" ~ "Cycle 12",
        AVISIT == "Week 24" ~ "Cycle 13",
        AVISIT == "Unscheduled 9.2" ~ "Cycle 15",
        AVISIT == "Unscheduled 9.3" ~ "Cycle 16",
        AVISIT == "Unscheduled 12.1" ~ "Cycle 17",
        AVISIT == "Unscheduled 13.1" ~ "Cycle 18",
        AVISIT == "Week 26" ~ "End Of Treatment",
        TRUE ~ as.character(AVISIT) # Other values remain unchanged
      )),
      AVISITN,
      .na_rm = FALSE
    ),

    # Demographic variables
    TRT01A = forcats::fct_reorder(TRT01A, TRT01AN, .na_rm = TRUE),
    TRTEMFL = as.factor(sample(c(NA, "Y"), dplyr::n(), replace = TRUE)),
    SEX = as.factor(dplyr::case_when(
      SEX == "F" ~ "Female",
      SEX == "M" ~ "Male"
    )),
    COUNTRY_DECODE = as.factor("United States of America"),
    RACE_DECODE = as.factor(dplyr::case_when(
      RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~
        "American Indian or Alaska Native",
      RACE == "ASIAN" ~ "Asian",
      RACE == "BLACK OR AFRICAN AMERICAN" ~ "Black or African American",
      RACE == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~
        "Native Hawaiian or other Pacific Islander",
      RACE == "WHITE" ~ "White",
      RACE == "MULTIPLE" ~ "Multiple",
      RACE == "NOT REPORTED" ~ "Not reported",
      RACE == "UNKNOWN" ~ "Unknown",
      RACE == "OTHER" ~ "Other"
    )),
    ETHNIC_DECODE = as.factor(dplyr::case_when(
      ETHNIC == "HISPANIC OR LATINO" ~ "Hispanic or Latino",
      ETHNIC == "NOT HISPANIC OR LATINO" ~ "Not Hispanic or Latino",
      ETHNIC == "NOT REPORTED" ~ "Not reported",
      ETHNIC == "UNKNOWN" ~ "Unknown"
    )),
    # Parameter coding
    PARAMCD = as.factor(case_when(
      PARAM == "Alkaline Phosphatase (U/L)" ~ "ALP",
      PARAM == "Potassium (mmol/L)" ~ "K",
      PARAM == "Cholesterol (mmol/L)" ~ "CHOL",
      PARAM == "Blood Urea Nitrogen (mmol/L)" ~ "LDL",
      PARAM == "Hematocrit (1)" ~ "HDL",
      PARAM == "Monocytes (10^9/L)" ~ "NEUT",
      PARAM == "Lymphocytes Abs (10^9/L)" ~ "TRIG",
      TRUE ~ as.character(PARAMCD)
    )),
    PARAM = as.factor(case_when(
      PARAMCD == "LDL" ~ "LDL Cholesterol (mmol/L)",
      PARAMCD == "HDL" ~ "HDL Cholesterol (mmol/L)",
      PARAMCD == "NEUT" ~ "Neutrophils (x10E9/L)",
      PARAMCD == "TRIG" ~ "Triglycerides (mmol/L)",
      PARAMCD == "WBC" ~ "Leukocytes (x10E9/L)",
      PARAMCD == "PLAT" ~ "Platelets (x10E9/L)",
      PARAMCD == "HGB" ~ "Hemoglobin (g/L)",
      TRUE ~ as.character(PARAM)
    )),
    AVAL = case_when(
      PARAMCD %in% c("LDL", "HDL") ~ runif(n(), min = 0.75, max = 2.5), # Generate random numbers between
      # 0.75 and 2.5 for LDL and HDL
      PARAMCD == "NEUT" ~ runif(n(), min = 0, max = 50), # Generate random numbers between 0 and 50 for NEUT
      PARAMCD == "TRIG" ~ runif(n(), min = 0.05, max = 5), # Generate random numbers between 0 and 50 for TRIG
      TRUE ~ AVAL # Keep the original AVAL for other cases
    ),
    PARCAT1 = as.factor(case_when(
      PARAM %in%
        c(
          "Alanine Aminotransferase (U/L)",
          "Albumin (g/L)",
          "Alkaline Phosphatase (U/L)",
          "Aspartate Aminotransferase (U/L)",
          "Bilirubin (umol/L)",
          "Calcium (mmol/L)",
          "Cholesterol (mmol/L)",
          "Corrected Calcium (mmol/L)",
          "Creatinine (umol/L)",
          "Direct Bilirubin (umol/L)",
          "HDL Cholesterol (mmol/L)",
          "Glucose (mmol/L)",
          "Indirect Bilirubin (umol/L)",
          "LDL Cholesterol (mmol/L)",
          "LDL Cholesterol (mmol/L) Calculated",
          "LDL Cholesterol (mmol/L) Direct",
          "Lactate Dehydrogenase (U/L)",
          "Potassium (mmol/L)",
          "Prostate Specific Antigen (ug/L)",
          "Protein (g/L)",
          "Serum Albumin (g/L)",
          "Sodium (mmol/L)",
          "Testosterone (nmol/L)",
          "Testosterone (nmol/L) Ultrasensitive Assay",
          "Thyrotropin (mIU/L)",
          "Thyroxine (nmol/L)",
          "Thyroxine, Free (pmol/L)",
          "Triglycerides (mmol/L)",
          "Triiodothyronine (nmol/L)"
        ) ~
        "CHEMISTRY",
      PARAM %in%
        c(
          "Blasts (x10E9/L)",
          "Hemoglobin (g/L)",
          "Leukocytes (x10E9/L)",
          "Neutrophils (x10E9/L)",
          "Neutrophils and Precursors (x10E9/L)",
          "Neutrophils, Segmented (x10E9/L)",
          "Platelets (x10E9/L)",
          "Prothrombin Intl. Normalized Ratio (RATIO)"
        ) ~
        "HEMATOLOGY",
      TRUE ~ NA_character_
    )),
    PARCAT2 = as.factor(case_when(
      PARAM == "Alanine Aminotransferase (U/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Albumin (g/L)" ~ "Test with FDA abnormality criteria defined",
      PARAM == "Alkaline Phosphatase (U/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Aspartate Aminotransferase (U/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Bilirubin (umol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Blasts (x10E9/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Calcium (mmol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Cholesterol (mmol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Creatinine (umol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Glucose (mmol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "HDL Cholesterol (mmol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "LDL Cholesterol (mmol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Triglycerides (mmol/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Hemoglobin (g/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Leukocytes (x10E9/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Neutrophils (x10E9/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "Platelets (x10E9/L)" ~
        "Test with FDA abnormality criteria defined",
      PARAM == "WBC Differential" ~
        "Test with FDA abnormality criteria defined",
      TRUE ~ NA_character_
    )),
    PARCAT3 = as.factor(case_when(
      PARAM == "Alanine Aminotransferase (U/L)" ~ "Liver biochemistry",
      PARAM == "Albumin (g/L)" ~ "Liver biochemistry",
      PARAM == "Alkaline Phosphatase (U/L)" ~ "Liver biochemistry",
      PARAM == "Aspartate Aminotransferase (U/L)" ~ "Liver biochemistry",
      PARAM == "Bilirubin (umol/L)" ~ "Liver biochemistry",
      PARAM == "Blasts (x10E9/L)" ~ "Liver biochemistry",
      PARAM == "Lactate Dehydrogenase (U/L)" ~ "Liver biochemistry",
      PARAM == "Protein (g/L)" ~ "Liver biochemistry",
      PARAM == "Calcium (mmol/L)" ~ "General chemistry",
      PARAM == "Creatinine (umol/L)" ~ "Kidney function",
      PARAM == "Creatinine Kinase (U/L)" ~ "Kidney function",
      PARAM == "Potassium (mmol/L)" ~ "General chemistry",
      PARAM == "Sodium (mmol/L)" ~ "General chemistry",
      PARAM == "LDL Cholesterol (mmol/L)" ~ "Lipids",
      PARAM == "Glucose (mmol/L)" ~ "General chemistry",
      PARAM == "HDL Cholesterol (mmol/L)" ~ "Lipids",
      PARAM == "Cholesterol (mmol/L)" ~ "Lipids",
      PARAM == "Triglycerides (mmol/L)" ~ "Lipids",
      PARAM == "Hemoglobin (g/L)" ~ "Complete blood count",
      PARAM == "WBC differential" ~ "Complete blood count",
      PARAM == "Platelets (x10E9/L)" ~ "Complete blood count",
      PARAM == "Leukocytes (x10E9/L)" ~ "Complete blood count",
      PARAM == "Neutrophils (x10E9/L)" ~ "WBC differential",
      PARAM == "Neutrophils, Segmented (x10E9/L)" ~ "WBC differential",
      PARAM == "Testosterone (nmol/L)" ~ "Endocrine",
      PARAM == "Thyroxine (nmol/L)" ~ "Endocrine",
      PARAM == "Thyrotropin (mIU/L)" ~ "Endocrine",
      PARAM == "Triiodothyronine (nmol/L)" ~ "Endocrine",
      TRUE ~ NA_character_
    )),
    PARCAT4 = as.factor(case_when(
      PARAM %in%
        c(
          "Alanine Aminotransferase (U/L)",
          "Alkaline Phosphatase (U/L)",
          "Aspartate Aminotransferase (U/L)",
          "Bilirubin (umol/L)",
          "Cholesterol (mmol/L)",
          "Creatinine (umol/L)",
          "Glucose (mmol/L)",
          "Hemoglobin (g/L)",
          "Leukocytes (x10E9/L)",
          "Neutrophils (x10E9/L)",
          "Neutrophils, Segmented (x10E9/L)",
          "Platelets (x10E9/L)",
          "Potassium (mmol/L)",
          "Protein (g/L)",
          "Prothrombin Intl. Normalized Ratio (RATIO)",
          "Sodium (mmol/L)",
          "Triglycerides (mmol/L)",
          "Neutrophils and Precursors (x10E9/L)",
          "Lactate Dehydrogenase (U/L)",
          "Direct Bilirubin (umol/L)",
          "HDL Cholesterol (mmol/L)"
        ) ~
        "Graded tests",
      TRUE ~ NA_character_
    )),
    PARCAT5 = as.factor(case_when(
      PARAM == "Alanine Aminotransferase (U/L)" ~ "Investigations",
      PARAM == "Albumin (g/L)" ~ "Metabolism and nutritional disorders",
      PARAM == "Alkaline Phosphatase (U/L)" ~ "Investigations",
      PARAM == "Aspartate Aminotransferase (U/L)" ~ "Investigations",
      PARAM == "Bilirubin (umol/L)" ~ "Investigations",
      PARAM == "Glucose (mmol/L)" ~ "Metabolism and nutritional disorders",
      PARAM == "Cholesterol (mmol/L)" ~ "Investigations",
      PARAM == "Creatinine (umol/L)" ~ "Investigations",
      PARAM == "HDL Cholesterol (mmol/L)" ~ "NA",
      PARAM == "Hemoglobin (g/L)" ~ "Investigations",
      PARAM == "Leukocytes (x10E9/L)" ~ "Investigations",
      PARAM == "Neutrophils (x10E9/L)" ~ "Investigations",
      PARAM == "Neutrophils, Segmented (x10E9/L)" ~ "Investigations",
      PARAM == "Platelets (x10E9/L)" ~ "Investigations",
      PARAM == "Potassium (mmol/L)" ~ "Metabolism and nutritional disorders",
      PARAM == "Protein (g/L)" ~ "Renal and urinary disorders",
      PARAM == "Prothrombin Intl. Normalized Ratio (RATIO)" ~ "Investigations",
      PARAM == "Sodium (mmol/L)" ~ "Metabolism and nutritional disorders",
      PARAM == "Triglycerides (mmol/L)" ~
        "Metabolism and nutritional disorders",
      TRUE ~ NA_character_
    )),
    PARCAT6 = as.factor(case_when(
      PARAM == "Direct Bilirubin (umol/L)" ~ "Bilirubin/Calcium Tests",
      PARAM == "Hemoglobin (g/L)" ~ "Blood and lymphatic system disorders",
      PARAM == "Leukocytes (x10E9/L)" ~ "Blood and lymphatic system disorders",
      PARAM == "Indirect Bilirubin (umol/L)" ~ "Bilirubin/Calcium Tests",
      PARAM == "Corrected Calcium (mmol/L)" ~ "Bilirubin/Calcium Tests",
      PARAM == "Blasts (x10E9/L)" ~ "Cellular Tests",
      PARAM == "Neutrophils and Precursors (x10E9/L)" ~ "Cellular Tests",
      TRUE ~ NA_character_
    )),

    # Sample variables
    MCRIT2ML = as.factor(sample(
      c(
        "Level 0",
        "Level 1 (>2.620 mmol/L)",
        "Level 1 (>5.5 mmol/L)",
        "Level 1 (>150 mmol/L)",
        "Level 2 (>13.0 10^9/L)",
        "Level 2 (>6.0 mmol/L)",
        "Level 3 (>6.5 mmol/L)",
        "Level 2 (>2.745 mmol/L)",
        "Level 1 (>10.8 10^9/L)",
        "Level 2 (fasting >=6.99 mmol/L or random >=11.10 mmol/L)",
        "Level 3 (>15.0 10^9/L)",
        NA
      ),
      size = n(),
      replace = TRUE
    )),
    MCRIT1ML = as.factor(sample(
      c(
        "Level 0",
        "Level 1 (125-135 g/L)",
        "Level 2 (>6.206 mmol/L)",
        "Level 1 (>5.172 mmol/L)",
        "Level 1 (<132 mmol/L)",
        "Level 1 (>1.694 mmol/L)",
        "Level 2 (>3.387 mmol/L)",
        "Level 1 (>1.5x ULN Enzyme U/L)",
        "Level 2 (<1.996 mmol/L)",
        NA
      ),
      size = n(),
      replace = TRUE
    )),
    MCRIT1MN = sample(c(0, 1, 2, 3, NaN), size = n(), replace = TRUE),
    MCRIT2MN = sample(c(0, 1, 2, 3, NaN), size = n(), replace = TRUE),
    # Multi-criteria variables
    MCRIT1 = case_when(
      PARAM == "Alanine Aminotransferase (U/L)" ~
        sample(c("Alanine Aminotransferase, high", NA), n(), replace = TRUE),
      PARAM == "Albumin (g/L)" ~
        sample(c("Albumin, low", NA), n(), replace = TRUE),
      PARAM == "Alkaline Phosphatase (U/L)" ~
        sample(c("Alkaline Phosphatase, high", NA), n(), replace = TRUE),
      PARAM == "Aspartate Aminotransferase (U/L)" ~
        sample(c("Aspartate Aminotransferase, high", NA), n(), replace = TRUE),
      PARAM == "Bilirubin (µmol/L)" ~
        sample(c("Bilirubin, high", NA), n(), replace = TRUE),
      PARAM == "Calcium (mmol/L)" ~
        sample(c("Calcium, low", NA), n(), replace = TRUE),
      PARAM == "Cholesterol (mmol/L)" ~
        sample(c("Cholesterol, high", NA), n(), replace = TRUE),
      PARAM == "Creatinine (µmol/L)" ~
        sample(c("Creatinine, low", NA), n(), replace = TRUE),
      PARAM == "Glucose (mmol/L)" ~
        sample(c("Glucose, low", NA), n(), replace = TRUE),
      PARAM == "HDL Cholesterol (mmol/L)" ~
        sample(c("HDL Cholesterol, males, low", NA), n(), replace = TRUE),
      PARAM == "Hemoglobin (g/L)" ~
        sample(c("Hemoglobin, male", NA), n(), replace = TRUE),
      PARAM == "LDL Cholesterol (mmol/L)" ~
        sample(c("LDL Cholesterol, high", NA), n(), replace = TRUE),
      PARAM == "Leukocytes (x10E9/L)" ~
        sample(c("Leukocytes, low", NA), n(), replace = TRUE),
      PARAM == "Neutrophils (x10E9/L)" ~
        sample(c("Neutrophils, low", NA), n(), replace = TRUE),
      PARAM == "Platelets (x10E9/L)" ~
        sample(c("Platelets, low", NA), n(), replace = TRUE),
      PARAM == "Potassium (mmol/L)" ~
        sample(c("Potassium, low", NA), n(), replace = TRUE),
      PARAM == "Protein (g/L)" ~
        sample(c("Protein, low", NA), n(), replace = TRUE),
      PARAM == "Sodium (mmol/L)" ~
        sample(c("Sodium, low", NA), n(), replace = TRUE),
      PARAM == "Triglycerides (mmol/L)" ~
        sample(c("Triglycerides, high"), n(), replace = TRUE),
      TRUE ~ NA_character_
    ),
    MCRIT2 = case_when(
      PARAM == "Calcium (mmol/L)" ~
        sample(c("Calcium, low", NA), n(), replace = TRUE),
      PARAM == "Glucose (mmol/L)" ~
        sample(c("Glucose, low", NA), n(), replace = TRUE),
      PARAM == "Leukocytes (x10E9/L)" ~
        sample(c("Leukocytes, low", NA), n(), replace = TRUE),
      PARAM == "Potassium (mmol/L)" ~
        sample(c("Potassium, low", NA), n(), replace = TRUE),
      PARAM == "Sodium (mmol/L)" ~
        sample(c("Sodium, low", NA), n(), replace = TRUE),
      TRUE ~ NA_character_ # Assign NA for other values not matching
    ),
    ATOXGR = case_when(
      ATOXGR == "0" ~ "0",
      ATOXGR == "1" ~ "1",
      ATOXGR == "-1" ~ "1",
      ATOXGR == "-2" ~ "2",
      ATOXGR == "-3" ~ "3",
      ATOXGR == "2" ~ "4",
      ATOXGR == "3" ~ "5",
      TRUE ~ NA_character_ # Other values can result in NA
    ),
    # Miscellaneous variables
    APOBLFL = as.factor(sample(c(NA, "Y"), dplyr::n(), replace = TRUE)),
    LBSTNRHQ = as.factor(sample(c(NA, "<"), dplyr::n(), replace = TRUE)),
    LBSTNRLQ = as.factor(sample(c(NA, "<"), dplyr::n(), replace = TRUE)),
    ATOXGRN = as.numeric(ATOXGR),
    ADTM = format(paste(ADT, "00:00"), format = "%Y-%m-%d %H:%M"),
    ATPT = strftime(ADTM, format = "%H:%M"),
    ATOXGRL = as.factor(sample(
      c(0, 1, 2, 3, 4, NaN),
      size = n(),
      replace = TRUE,
      prob = c(0.618, 0.2, 0.1, 0.05, 0.005, 0.067)
    )),
    TR01SDT = sample(
      seq(
        min(as.Date(TRTSDT), na.rm = TRUE),
        max(as.Date(TRTSDT), na.rm = TRUE),
        by = "day"
      ),
      length(TRTEDT),
      replace = TRUE
    ),
    TR01EDT = sample(
      seq(
        min(as.Date(TRTEDT), na.rm = TRUE),
        max(as.Date(TRTEDT), na.rm = TRUE),
        by = "day"
      ),
      length(TRTEDT),
      replace = TRUE
    )
  )

  # Apply admiral::restrict_derivation for ANL03FL
  gen <- admiral::restrict_derivation(
    gen,
    derivation = admiral::derive_var_extreme_flag,
    args = admiral::params(
      by_vars = rlang::syms(c("USUBJID", "PARAMCD", "AVISIT")),
      order = rlang::syms(c("AVAL", "ADT", "ADY")),
      new_var = ANL03FL,
      true_value = "Y",
      false_value = NA,
      mode = "last"
    ),
    filter = AVISIT != "Screening"
  )

  # Apply admiral::restrict_derivation for ANL04FL
  gen <- admiral::restrict_derivation(
    gen,
    derivation = admiral::derive_var_extreme_flag,
    args = admiral::params(
      by_vars = rlang::syms(c("USUBJID", "PARAMCD", "AVISIT")),
      order = rlang::syms(c("AVAL", "ADT", "ADY")),
      new_var = ANL04FL,
      true_value = "Y",
      false_value = NA,
      mode = "last"
    ),
    filter = AVISIT != "Screening"
  )

  # Apply admiral::restrict_derivation for ANL05FL
  gen <- admiral::restrict_derivation(
    gen,
    derivation = admiral::derive_var_extreme_flag,
    args = admiral::params(
      by_vars = rlang::syms(c("USUBJID", "PARAMCD", "AVISIT")),
      order = rlang::syms(c("AVAL", "ADT", "ADY")),
      new_var = ANL05FL,
      true_value = "Y",
      false_value = NA,
      mode = "last"
    ),
    filter = AVISIT != "Screening"
  ) %>%
    mutate(
      ANL06FL = ANL05FL,
      ANL07FL = ANL05FL,
      ANL08FL = ANL05FL,
      ANL09FL = ANL05FL,
      ANL10FL = ANL05FL,
      ANL14FL = ANL05FL,
      ANL15FL = ANL05FL,
      ANL16FL = ANL05FL
    )

  # Additional labels for new variables not in the source dataset
  additional_labels <- list(
    ANL02FL = "Analysis Record Flag 02-Analysis Value",
    ANL03FL = "Analysis Record Flag 03 - Protocol Visit",
    ANL04FL = "Analysis Flag 04",
    ANL05FL = "Analysis Flag 05",
    ANL06FL = "Analysis Flag 06",
    ANL07FL = "Analysis Flag 07",
    ANL08FL = "Analysis Flag 08",
    ANL09FL = "Analysis Flag 09",
    ANL10FL = "Analysis Flag 10",
    ANL14FL = "Analysis Flag 14",
    ANL15FL = "Analysis Flag 15",
    ANL16FL = "Analysis Flag 16",
    APOBLFL = "Post-Baseline Record Flag",
    LBSTNRHQ = "Reference Limit Higher",
    LBSTNRLQ = "Reference Limit Lower",
    TRT01P = "Planned Treatment for Period 01",
    TRT01A = "Actual Treatment for Period 01",
    ADT = "Analysis Date",
    AVAL = "Analysis Value",
    BASE = "Baseline Value",
    CHG = "Change from Baseline",
    PCHG = "Percent Change from Baseline",
    SEX = "Sex",
    RACE = "Race",
    COUNTRY_DECODE = "Country",
    PARCAT1 = "Parameter Category 1",
    PARCAT2 = "Parameter Category 2",
    PARCAT3 = "Parameter Category 3",
    PARCAT4 = "Parameter Category 4",
    PARCAT5 = "Parameter Category 5",
    PARCAT6 = "Parameter Category 6",
    ATOXGR = "Analysis Toxicity Grade",
    MCRIT1 = "Analysis Multi-Response Criterion 1",
    MCRIT2 = "Analysis Multi-Response Criterion 2",
    MCRIT1ML = "Multi-Response Criterion 1 Evaluation",
    MCRIT2ML = "Multi-Response Criterion 2 Evaluation",
    MCRIT1MN = "Multi-Response Criterion 1 Eval (N)",
    MCRIT2MN = "Multi-Response Criterion 2 Eval (N)",
    TRTSDT = "Date of First Exposure to Treatment",
    TRTEDT = "Date of Last Exposure to Treatment",
    AVISIT = "Analysis Visit",
    AVISITN = "Analysis Visit (N)",
    TRT01PN = "Planned Treatment for Period 01 (N)",
    TRT01AN = "Actual Treatment for Period 01 (N)",
    ADYM = "Days from Treatment Start",
    ATPT = "Analysis Timepoint",
    AVALC = "Analysis Value (C)",
    AVALU = "Analysis Value - Units",
    TRTEMFL = "Treatment Emergent Analysis Flag",
    RACE_DECODE = "Race Description",
    ETHNIC_DECODE = "Ethnicity Description",
    ATOXGRN = "Analysis Toxicity Grade (Numeric)",
    ATOXGRL = "Analysis Toxicity Grade Low",
    ADTM = "Analysis Date/Time",
    TRT01SDT = "Start Date of Planned Treatment for Period 01",
    TRT01EDT = "End Date of Planned Treatment for Period 01",
    TR01SDT = "Start Date of Treatment for Period 01",
    TR01EDT = "End Date of Treatment for Period 01"
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
adlb <- gen_adlb()
