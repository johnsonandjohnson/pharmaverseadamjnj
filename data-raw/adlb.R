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

  gen$TRT01P <- as.factor(gen$TRT01P)
  gen$TRT01A <- as.factor(gen$TRT01A)

  # Added Toxicity grade lookup table as per latest version
  toxterm_lookup <- tibble::tribble(
    ~PARAMCD, ~ATOXDSCL, ~ATOXDSCH, ~ATOXDIR,
    "ALB", "Albumin, low", NA, "LOW",
    "ALP", NA, "Alkaline Phosphatase, high", "HIGH",
    "ALT", NA, "Alanine Aminotransferase, high", "HIGH",
    "AMYLASE", NA, "Amylase, high", "HIGH",
    "APTT", NA, "Activated Partial Thromboplastin Time, high", "HIGH",
    "AST", NA, "Aspartate Aminotransferase, high", "HIGH",
    "BILI", NA, "Bilirubin, high", "HIGH",
    "CACRALB", NA, "Calcium Corrected, high", "HIGH",
    "CACRALB", "Calcium Corrected, low", NA, "LOW",
    "CACR", NA, "Calcium Corrected, high", "HIGH",
    "CACR", "Calcium Corrected, low", NA, "LOW",
    "CAION", NA, "Calcium, Ionized, high", "HIGH",
    "CAION", "Calcium, Ionized, low", NA, "LOW",
    "CD4", "CD4, low", NA, "LOW",
    "CHOL", NA, "Cholesterol, high", "HIGH",
    "CK", NA, "Creatine Kinase, high", "HIGH",
    "CREAT", NA, "Creatinine, high", "HIGH",
    "FIBRINO", "Fibrinogen, decreased", NA, "LOW",
    "GGT", NA, "Gamma Glutamyl Transferase, high", "HIGH",
    "GLUC", "Glucose, low", NA, "LOW",
    "HAPTOG", "Haptoglobin, low", NA, "LOW",
    "HGB", NA, "Hemoglobin, high", "HIGH",
    "HGB", "Hemoglobin, low", NA, "LOW",
    "INR", NA, "Prothrombin Intl. Normalized Ratio, high", "HIGH",
    "K", NA, "Potassium, high", "HIGH",
    "K", "Potassium, low", NA, "LOW",
    "LIPASET", NA, "Lipase, high", "HIGH",
    "LYM", NA, "Lymphocytes, high", "HIGH",
    "LYM", "Lymphocytes, low", NA, "LOW",
    "MG", NA, "Magnesium, high", "HIGH",
    "MG", "Magnesium, low", NA, "LOW",
    "NEUT", "Neutrophils, low", NA, "LOW",
    "NEUTSG", "Neutrophils, low", NA, "LOW",
    "NEUTSGB", "Neutrophils, low", NA, "LOW",
    "NEUTSGBP", "Neutrophils, low", NA, "LOW",
    "PH", NA, "pH, high", "HIGH",
    "PH", "pH, low", NA, "LOW",
    "PLAT", "Platelets, low", NA, "LOW",
    "PROT", NA, "Urinary Protein, high", "HIGH",
    "PROTCRT", NA, "Urinary Protein, high", "HIGH",
    "SODIUM", NA, "Sodium, high", "HIGH",
    "SODIUM", "Sodium, low", NA, "LOW",
    "TRIG", NA, "Triglycerides, high", "HIGH",
    "WBC", NA, "Leukocytes, high", "HIGH",
    "WBC", "Leukocytes, low", NA, "LOW"
  )

  # Adding latest ATOXDSCL & ATOXDSCH based on LBTESTCD from toxterm lookup table and removed the existing variables
  gen <- gen |>
    select(-c(ATOXDSCL, ATOXDSCH)) |>
    admiral::derive_vars_merged(
      dataset_add = dplyr::filter(toxterm_lookup, ATOXDIR == "LOW") |>
        dplyr::select(PARAMCD, ATOXDSCL),
      by_vars = exprs(PARAMCD)
    ) |>
    admiral::derive_vars_merged(
      dataset_add = dplyr::filter(toxterm_lookup, ATOXDIR == "HIGH") |>
        dplyr::select(PARAMCD, ATOXDSCH),
      by_vars = exprs(PARAMCD)
    )

  gen <- dplyr::mutate(
    gen,
    # adjust AVAL and ANRHI for ALKPH to ensure ratio > 3 when calculated on-the-fly
    AVAL = ifelse(
      PARAMCD == "ALKPH" & USUBJID == unique(USUBJID)[1] & ADY > 0,
      30, # Use a high AVAL
      AVAL
    ),
    ANRHI = ifelse(
      PARAMCD == "ALKPH" & USUBJID == unique(USUBJID)[1] & ADY > 0,
      5, # ensure AVAL/ANRHI = 30/5 = 6 > 3
      ANRHI
    ),
    # Treatment and arm variables
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
    TRT01AN = dplyr::case_when(
      TRT01A == "Xanomeline High Dose" ~ 1,
      TRT01A == "Xanomeline Low Dose" ~ 2,
      TRT01A == "Placebo" ~ 3
    ),

    # Analysis values
    AVAL = round(AVAL, 4),
    AVALC = as.character(AVAL),
    AVALU = LBSTRESU,
    ANL02FL = "Y",
    AVISIT = case_when(
      ABLFL == "Y" ~ "Baseline",
      toupper(AVISIT) == "BASELINE" & (is.na(ABLFL) | ABLFL != "Y") ~ VISIT,
      TRUE ~ AVISIT
    ),
    # Visit variables
    AVISITN = case_when(
      toupper(AVISIT) == "BASELINE" ~ 1,
      toupper(AVISIT) == "WEEK 2" ~ 2,
      toupper(AVISIT) == "WEEK 4" ~ 3,
      toupper(AVISIT) == "WEEK 6" ~ 4,
      toupper(AVISIT) == "WEEK 8" ~ 5,
      toupper(AVISIT) == "WEEK 12" ~ 6,
      toupper(AVISIT) == "WEEK 16" ~ 7,
      toupper(AVISIT) == "WEEK 20" ~ 8,
      toupper(AVISIT) == "WEEK 24" ~ 9,
      toupper(AVISIT) == "WEEK 26" ~ 10
    ),
    AVISIT = fct_reorder(
      as.factor(case_when(
        toupper(AVISIT) == "BASELINE" ~ "Baseline",
        toupper(AVISIT) == "WEEK 2" ~ "Cycle 02",
        toupper(AVISIT) == "WEEK 4" ~ "Cycle 03",
        toupper(AVISIT) == "WEEK 6" ~ "Cycle 04",
        toupper(AVISIT) == "WEEK 8" ~ "Cycle 05",
        toupper(AVISIT) == "WEEK 12" ~ "Cycle 12",
        toupper(AVISIT) == "WEEK 16" ~ "Cycle 16",
        toupper(AVISIT) == "WEEK 20" ~ "Cycle 20",
        toupper(AVISIT) == "WEEK 24" ~ "Cycle 24",
        toupper(AVISIT) == "WEEK 26" ~ "End Of Treatment",
        TRUE ~ stringr::str_to_sentence(as.character(AVISIT))
      )),
      AVISITN,
      .na_rm = FALSE
    ),

    # Demographic variables
    TRT01A = forcats::fct_reorder(TRT01A, TRT01AN, .na_rm = TRUE),
    TRTEMFL = as.factor(sample(c(NA, "Y"), dplyr::n(), replace = TRUE)),
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

    # MCRIT1ML: criterion 1 evaluation level based on PARAM
    MCRIT1ML = as.factor(case_when(
      PARAM == "Alanine Aminotransferase (U/L)" ~
        sample(c(
          "Level 0", "Level 1 (>1.5x ULN U/L)",
          "Level 2 (>3.0x ULN U/L)", "Level 3 (>5.0x ULN U/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Albumin (g/L)" ~
        sample(c("Level 0", "Level 1 (25-35 g/L)", "Level 2 (<25 g/L)", NA), n(), replace = TRUE),
      PARAM == "Alkaline Phosphatase (U/L)" ~
        sample(c(
          "Level 0", "Level 1 (>1.5x ULN U/L)",
          "Level 2 (>3.0x ULN U/L)", "Level 3 (>5.0x ULN U/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Aspartate Aminotransferase (U/L)" ~
        sample(c(
          "Level 0", "Level 1 (>1.5x ULN U/L)",
          "Level 2 (>3.0x ULN U/L)", "Level 3 (>5.0x ULN U/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Bilirubin (umol/L)" ~
        sample(c("Level 0", "Level 1 (>1.5x ULN umol/L)", "Level 2 (>3.0x ULN umol/L)", NA), n(), replace = TRUE),
      PARAM == "Calcium (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (<2.0 mmol/L)",
          "Level 2 (<1.75 mmol/L)", "Level 3 (<1.5 mmol/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Cholesterol (mmol/L)" ~
        sample(c("Level 0", "Level 1 (>5.172 mmol/L)", "Level 2 (>6.206 mmol/L)", NA), n(), replace = TRUE),
      PARAM == "Creatinine (umol/L)" ~
        sample(c("Level 0", "Level 1 (>1.5x ULN umol/L)", "Level 2 (>3.0x ULN umol/L)", NA), n(), replace = TRUE),
      PARAM == "Glucose (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (<3.0 mmol/L)",
          "Level 2 (fasting >=6.99 mmol/L or random >=11.10 mmol/L)", NA
        ), n(), replace = TRUE),
      PARAM == "HDL Cholesterol (mmol/L)" ~
        sample(c("Level 0", "Level 1 (<1.036 mmol/L)", "Level 2 (<0.777 mmol/L)", NA), n(), replace = TRUE),
      PARAM == "Hemoglobin (g/L)" ~
        sample(c(
          "Level 0", "Level 1 (100-119 g/L)",
          "Level 2 (80-99 g/L)", "Level 3 (<80 g/L)", NA
        ), n(), replace = TRUE),
      PARAM == "LDL Cholesterol (mmol/L)" ~
        sample(c("Level 0", "Level 1 (>3.367 mmol/L)", "Level 2 (>4.137 mmol/L)", NA), n(), replace = TRUE),
      PARAM == "Leukocytes (x10E9/L)" ~
        sample(c(
          "Level 0", "Level 1 (<3.0 10^9/L)",
          "Level 2 (<2.0 10^9/L)", "Level 3 (<1.0 10^9/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Neutrophils (x10E9/L)" ~
        sample(c(
          "Level 0", "Level 1 (<1.5 10^9/L)",
          "Level 2 (<1.0 10^9/L)", "Level 3 (<0.5 10^9/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Platelets (x10E9/L)" ~
        sample(c(
          "Level 0", "Level 1 (<LLN-75.0 10^9/L)",
          "Level 2 (<75.0-50.0 10^9/L)", "Level 3 (<50.0 10^9/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Potassium (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (<3.0 mmol/L)",
          "Level 2 (<2.5 mmol/L)", "Level 3 (<2.0 mmol/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Protein (g/L)" ~
        sample(c("Level 0", "Level 1 (1+ g/L)", "Level 2 (2+ g/L)", "Level 3 (3+ g/L)", NA), n(), replace = TRUE),
      PARAM == "Sodium (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (<130 mmol/L)",
          "Level 2 (<125 mmol/L)", "Level 3 (<120 mmol/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Triglycerides (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (>1.694 mmol/L)",
          "Level 2 (>2.260 mmol/L)", "Level 3 (>5.650 mmol/L)", NA
        ), n(), replace = TRUE),
      TRUE ~ NA_character_
    )),
    # MCRIT2ML: criterion 2 evaluation level based on PARAM (only PARAMs with MCRIT2 defined)
    MCRIT2ML = as.factor(case_when(
      PARAM == "Calcium (mmol/L)" ~
        sample(c("Level 0", "Level 1 (>2.620 mmol/L)", "Level 2 (>2.745 mmol/L)", NA), n(), replace = TRUE),
      PARAM == "Glucose (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (>5.5 mmol/L)",
          "Level 2 (>6.0 mmol/L)", "Level 3 (>6.5 mmol/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Leukocytes (x10E9/L)" ~
        sample(c(
          "Level 0", "Level 1 (>10.8 10^9/L)",
          "Level 2 (>13.0 10^9/L)", "Level 3 (>15.0 10^9/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Potassium (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (>5.0 mmol/L)",
          "Level 2 (>5.5 mmol/L)", "Level 3 (>6.0 mmol/L)", NA
        ), n(), replace = TRUE),
      PARAM == "Sodium (mmol/L)" ~
        sample(c(
          "Level 0", "Level 1 (>145 mmol/L)",
          "Level 2 (>150 mmol/L)", "Level 3 (>155 mmol/L)", NA
        ), n(), replace = TRUE),
      TRUE ~ NA_character_
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
      TRUE ~ NA_character_
    ),
    # NOTE: CRIT1, CRIT2, CRIT1FL, CRIT2FL are temporarily derived from MCRITy/MCRITyML
    # as the markedly abnormal criteria definitions file is currently missing.
    # These should be updated to reflect the actual criterion definitions and flags
    # from the markedly abnormal file once it is available.
    CRIT1 = MCRIT1,
    CRIT2 = MCRIT2,
    CRIT1FL = MCRIT1ML,
    CRIT2FL = MCRIT2ML,
    ATOXGR = as.character(sample(
      c("0", "1", "2", "3", NA_character_),
      size = n(),
      replace = TRUE,
      prob = c(0.5, 0.25, 0.15, 0.07, 0.03)
    )),
    # Miscellaneous variables
    APOBLFL = as.factor(dplyr::if_else(
      (is.na(ABLFL) | ABLFL != "Y") & !is.na(ADT) & !is.na(TRTSDT) & as.Date(ADT) >= as.Date(TRTSDT),
      "Y",
      NA_character_
    )),
    LBSTNRHQ = as.factor(sample(c(NA, "<"), dplyr::n(), replace = TRUE)),
    LBSTNRLQ = as.factor(sample(c(NA, "<"), dplyr::n(), replace = TRUE)),
    ATOXGRN = as.numeric(ATOXGR),
    ADTM = format(paste(ADT, "00:00"), format = "%Y-%m-%d %H:%M"),
    ATPT = "BEFORE TREATMENT",
    ATOXGRL = as.factor(dplyr::if_else(!is.na(ATOXDSCL), as.numeric(ATOXGR), NA_real_)),
    ATOXGRH = as.factor(dplyr::if_else(!is.na(ATOXDSCH), as.numeric(ATOXGR), NA_real_)),
    # Add LBCLSIG variable with values "N" and "Y"
    LBCLSIG = as.factor(sample(c("N", "Y"), size = n(), replace = TRUE, prob = c(0.7, 0.3))),
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
    ),
    LBSPEC = dplyr::case_when(
      PARAMCD == "GLUC" ~ "PLASMA"
    ),
    LBFAST = dplyr::case_when(
      PARAMCD == "GLUC" ~ "Y"
    ),
    LBNAM = sample(c("CENTRAL", "LOCAL"), n(), replace = TRUE, prob = c(0.85, 0.15))
  )

  # Baseline Toxicity derivation
  gen <- gen |>
    select(-c(BTOXGRH, BTOXGRL, BTOXGR)) |>
    derive_var_base(
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      source_var = ATOXGRL,
      new_var = BTOXGRL
    ) |>
    derive_var_base(
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      source_var = ATOXGRH,
      new_var = BTOXGRH
    ) |>
    derive_var_base(
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      source_var = ATOXGR,
      new_var = BTOXGR
    )

  # ANL01FL should be flaged for non missing records that is AVAL or AVALC not missing
  # Apply admiral::restrict_derivation for ANL01FL
  gen <- gen |>
    select(-(ANL01FL)) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(USUBJID, PARAMCD, AVISIT),
        order = exprs(ADT, AVAL),
        new_var = ANL01FL,
        mode = "last"
      ),
      filter = !is.na(AVAL) & !is.na(AVALC)
    )

  # Derivation of ANL02FL flag for scheduled visits
  gen <- gen |>
    mutate(
      ANL02FL = dplyr::case_when(
                                 (grepl("Cycle", AVISIT) |
                                    grepl("End Of Treatment", AVISIT) |
                                       grepl("Baseline", AVISIT)) &
                                 ANL01FL == "Y" & is.na(DTYPE) ~ "Y",
                                 TRUE ~ NA_character_)
    )

  # Apply admiral::restrict_derivation for ANL03FL
  gen <- admiral::restrict_derivation(
    gen,
    derivation = admiral::derive_var_extreme_flag,
    args = admiral::params(
      by_vars = rlang::syms(c("USUBJID", "PARAMCD")),
      order = rlang::syms(c("AVAL", "AVISIT", "ADT", "ADY")),
      new_var = ANL03FL,
      true_value = "Y",
      false_value = NA,
      mode = "last"
    ),
    filter = ONTRTFL == "Y" & ANL01FL == "Y"
  )

  # Apply admiral::restrict_derivation for ANL04FL
  gen <- admiral::restrict_derivation(
    gen,
    derivation = admiral::derive_var_extreme_flag,
    args = admiral::params(
      by_vars = rlang::syms(c("USUBJID", "PARAMCD")),
      order = rlang::syms(c("ATOXGRL", "AVISIT", "ADT", "ADY")),
      new_var = ANL04FL,
      true_value = "Y",
      false_value = NA,
      mode = "first"
    ),
    filter = ONTRTFL == "Y" & ANL01FL == "Y"
  )

  # Apply admiral::restrict_derivation for ANL05FL
  gen <- admiral::restrict_derivation(
    gen,
    derivation = admiral::derive_var_extreme_flag,
    args = admiral::params(
      by_vars = rlang::syms(c("USUBJID", "PARAMCD")),
      order = rlang::syms(c("ATOXGRL", "AVISIT", "ADT", "ADY")),
      new_var = ANL05FL,
      true_value = "Y",
      false_value = NA,
      mode = "last"
    ),
    filter = ONTRTFL == "Y" & ANL01FL == "Y"
  ) |>
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

  # Replace x10E with 10^ in PARAM labels
  gen <- dplyr::mutate(gen, PARAM = as.factor(gsub("x10E", "10^", as.character(PARAM))))

  source(file.path("data-raw", "adsl.R"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), "USUBJID")

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c(
    "TRT01A",
    "SAFFL",
    "STUDYID",
    "AGE",
    "SEX",
    "RACE"
  )

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl |>
    dplyr::select(USUBJID, dplyr::all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID")

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
    LBCLSIG = "Clinically Significant",
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
    PARCAT1 = "Parameter Category 1",
    PARCAT2 = "Parameter Category 2",
    PARCAT3 = "Parameter Category 3",
    PARCAT4 = "Parameter Category 4",
    PARCAT5 = "Parameter Category 5",
    PARCAT6 = "Parameter Category 6",
    ATOXGR = "Analysis Toxicity Grade",
    MCRIT1 = "Analysis Multi-Response Criterion 1",
    MCRIT2 = "Analysis Multi-Response Criterion 2",
    CRIT1 = "Analysis Criterion 1",
    CRIT2 = "Analysis Criterion 2",
    CRIT1FL = "Criterion 1 Evaluation Flag",
    CRIT2FL = "Criterion 2 Evaluation Flag",
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
    ATOXGRN = "Analysis Toxicity Grade (Numeric)",
    ATOXGRL = "Analysis Toxicity Grade Low",
    ATOXGRH = "Analysis Toxicity Grade High",
    ADTM = "Analysis Date/Time",
    TRT01SDT = "Start Date of Planned Treatment for Period 01",
    TRT01EDT = "End Date of Planned Treatment for Period 01",
    TR01SDT = "Start Date of Treatment for Period 01",
    TR01EDT = "End Date of Treatment for Period 01",
    LBSPEC = "Specimen Type",
    LBFAST = "Fasting Status",
    LBNAM = "Laboratory Name"
  )

  # Handle NA values and convert characters to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Restore labels
  gen <- restore_labels(
    df = gen,
    orig_df = raw,
    additional_labels = additional_labels
  )

  # Ensure uniqueness of records per subject/parameter/visit
  gen <- gen |>
    group_by(USUBJID, PARAMCD, AVISIT) |>
    slice(1) |>
    ungroup()

  return(gen)
}

# Generate the dataset
adlb <- gen_adlb()
