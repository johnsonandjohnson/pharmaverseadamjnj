# Generate ADCM dataset

# Load necessary libraries
library(dplyr)
library(pharmaverseadam)
library(formatters)

# Source utility functions
source(file.path("data-raw", "helpers.R"))


# Generate ADCM dataset
gen_adcm <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Get source data
  raw <- pharmaverseadam::adcm

  gen <- raw
  gen$CMLVL1 <- as.factor(gen$CMCLAS)
  gen$CMLVL2 <- as.factor(gen$CMCLAS)
  gen$CMLVL3 <- as.factor(gen$CMCLAS)
  gen$CMLVL4 <- as.factor(gen$CMCLAS)
  gen$CMBASPRF <- as.factor(gen$CMDECOD)
  gen$CMPRESP <- as.factor("Y")
  gen$CMOCCUR <- as.factor("Y")
  gen$CMINDCSP <- as.factor(gen$CMINDC)
  gen$CMDOSTXT <- as.factor(gen$CMDOSE)
  gen$CMENRF <- as.factor(dplyr::case_when(
    !is.na(gen$CMENRTPT) ~ "AFTER",
    .default = NA
  ))
  gen$PREFL <- as.factor(dplyr::case_when(
    gen$ASTDT < gen$TRTSDT ~ "Y",
    .default = NA
  ))
  gen$ONTRTFL <- as.factor(dplyr::case_when(
    gen$TRTSDT < gen$ASTDT & gen$ASTDT < gen$TRTEDT + 30 ~ "Y",
    .default = NA
  ))
  gen$FUPFL <- as.factor(dplyr::case_when(
    gen$ASTDT > gen$TRTEDT + 30 ~ "Y",
    .default = NA
  ))
  gen$CQ01NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "HYDROCORTISONE" ~ "Steroid",
    .default = NA
  ))
  gen$CQ02NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "ACETYLSALICYLIC ACID" ~ "Acetylsalicylic Acid",
    .default = NA
  ))
  gen$CQ03NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "AMLODIPINE" ~ "Amlodipine",
    .default = NA
  ))
  gen$CQ04NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "NIZATIDINE" ~ "Nizatidine",
    .default = NA
  ))
  gen$CQ05NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "NIFEDIPINE" ~ "Nifedipine",
    .default = NA
  ))
  gen$CQ06NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "FUROSEMIDE" ~ "Furosemide",
    .default = NA
  ))
  gen$CQ07NAM <- as.factor(dplyr::case_when(
    gen$CMDECOD == "SALBUTAMOL SULFATE" ~ "Salbutamol Sulfate",
    .default = NA
  ))

  # remove orphan categories
  gen <- filter(
    gen,
    !CMBASPRF %in%
      c("AMLODIPINE", "NIZATIDINE", "FUROSEMIDE", "SALBUTAMOL SULFATE")
  )

  # Load the adsl dataset
  source(file.path("data-raw", "adsl.R"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), "USUBJID")

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c("TRT01A", "SAFFL", "TRTSDT")

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl %>%
    select(USUBJID, all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    message("Dropping shared vars from raw: ", paste(shared, collapse = ", "))
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID")

  # Additional labels for new variables not in the source dataset
  additional_labels <- list(
    # Concomitant medication variables
    CMLVL1 = "Preferred ATC Text for ATC Level 1",
    CMLVL2 = "Preferred ATC Text for ATC Level 2",
    CMLVL3 = "Preferred ATC Text for ATC Level 3",
    CMLVL4 = "Preferred ATC Text for ATC Level 4",
    CMBASPRF = "Base Preferred Term",
    CMPRESP = "CM Pre-specified",
    CMOCCUR = "CM Occurrence",
    CMINDCSP = "Indication Specification",
    CMDOSTXT = "Dose Description",
    CMENRTPT = "End Relative to Reference Time Point",
    CMENRF = "End Relative to Reference Period",

    # Custom query variables
    CQ01NAM = "Customized Query 01 Name",
    CQ02NAM = "Customized Query 02 Name",
    CQ03NAM = "Customized Query 03 Name",
    CQ04NAM = "Customized Query 04 Name",
    CQ05NAM = "Customized Query 05 Name",
    CQ06NAM = "Customized Query 06 Name",
    CQ07NAM = "Customized Query 07 Name"
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


adcm <- gen_adcm()
