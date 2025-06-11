# Generate ADTTESAF dataset

# Load necessary libraries
library(dplyr)
library(forcats)
library(pharmaverseadam)
library(formatters)
library(labelled)

# Source utility functions
source(file.path("data-raw", "helpers.R"))

# Generate ADTTESAF dataset
gen_adttesaf <- function(seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Start with pharmaverseadam::adtte_onco for structure
  raw <- pharmaverseadam::adtte_onco

  # Get unique subject IDs from the base data
  unique_subjects <- unique(raw$USUBJID)

  # Define the specific parameter values needed
  param_mapping <- data.frame(
    PARAMCD = c(
      "TTDCEOSM",
      "TTLSTFUD",
      "TTLSTFUM",
      "TTLSTFUW",
      "TTDCEOSD",
      "TTAELPTD",
      "TTDCEOSW",
      "TTDCEOTD"
    ),
    PARAM = c(
      "Time to Discontinuation from Study (months)",
      "Time To Last Follow-Up (days)",
      "Time To Last Follow-Up (months)",
      "Time To Last Follow-Up (weeks)",
      "Time to Discontinuation from Study (days)",
      "Time To Adverse Events Leading to Permanent Treatment Discontinuation (days)",
      "Time to Discontinuation from Study (weeks)",
      "Time to Discontinuation from Treatment (days)"
    ),
    stringsAsFactors = FALSE
  )

  # Create empty list to store records
  records <- list()
  record_index <- 1

  # For each subject, create one record for each parameter
  for (subj in unique_subjects) {
    for (i in seq_len(nrow(param_mapping))) {
      # Use the base structure but replace key fields
      record <- raw[1, ] # Use first row as template
      record$USUBJID <- subj
      record$PARAMCD <- param_mapping$PARAMCD[i]
      record$PARAM <- param_mapping$PARAM[i]

      # Generate appropriate AVAL based on parameter
      if (grepl("month", param_mapping$PARAM[i], ignore.case = TRUE)) {
        record$AVAL <- round(runif(1, 0, 36), 1) # 0-36 months
      } else if (grepl("week", param_mapping$PARAM[i], ignore.case = TRUE)) {
        record$AVAL <- round(runif(1, 0, 156), 0) # 0-156 weeks (3 years)
      } else {
        # days
        record$AVAL <- round(runif(1, 0, 1095), 0) # 0-1095 days (3 years)
      }

      # Set CNSR with approximately 70% censoring rate
      record$CNSR <- sample(c(0, 1), 1, prob = c(0.3, 0.7))

      # Add to records list
      records[[record_index]] <- record
      record_index <- record_index + 1
    }
  }

  # Combine all records
  gen <- dplyr::bind_rows(records)

  # Select only needed columns from base data
  gen <- dplyr::select(
    gen,
    USUBJID,
    PARAMCD,
    PARAM,
    AVAL,
    CNSR
  )

  # Add additional
  gen$STARTDT <- ifelse(
    runif(nrow(gen)) < 0.1,
    NA,
    as.Date("2016-01-01") + sample(0:730, nrow(gen), replace = TRUE)
  )
  gen$ADT <- ifelse(
    is.na(gen$STARTDT),
    NA,
    gen$STARTDT + sample(0:1000, nrow(gen), replace = TRUE)
  )

  source(file.path("data-raw", "adsl.R"))

  # Drop any variables shared by gen and ADSL (except key)
  shared <- setdiff(intersect(names(gen), names(adsl)), "USUBJID")

  # Variables to keep exclusively from ADSL
  to_keep_from_adsl <- c("TRT01A", "SAFFL")

  # Select only the key and the 'to_keep' variables from ADSL
  adsl_subset <- adsl %>%
    select(USUBJID, all_of(to_keep_from_adsl))

  if (length(shared) > 0) {
    message("Dropping shared vars from raw: ", paste(shared, collapse = ", "))
    gen <- dplyr::select(gen, -dplyr::any_of(shared))
  }

  gen <- dplyr::left_join(gen, adsl_subset, by = "USUBJID")

  # Convert character vars to factors
  gen <- df_na(gen, char_as_factor = TRUE)

  # Define labels for key variables
  additional_labels <- list(
    AVAL = "Analysis Value",
    CNSR = "Censor",
    USUBJID = "Unique Subject Identifier",
    PARAMCD = "Parameter Code",
    PARAM = "Parameter",
    FASFL = "Full Analysis Set Flag",
    TTDCEOSM = "Time to Discontinuation from Study (months)",
    TTLSTFUD = "Time To Last Follow-Up (days)",
    TTLSTFUM = "Time To Last Follow-Up (months)",
    TTLSTFUW = "Time To Last Follow-Up (weeks)",
    TTDCEOSD = "Time to Discontinuation from Study (days)",
    TTAELPTD = "Time To Adverse Events Leading to Permanent Treatment Discontinuation (days)",
    TTDCEOSW = "Time to Discontinuation from Study (weeks)",
    TTDCEOTD = "Time to Discontinuation from Treatment (days)",
    STARTDT = "Start Date",
    ADT = "Analysis Date"
  )

  gen <- restore_labels(
    df = gen,
    orig_df = raw,
    additional_labels = additional_labels
  )

  return(gen)
}

# Generate the dataset
adttesaf <- gen_adttesaf()
