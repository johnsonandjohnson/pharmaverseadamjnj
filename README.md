# pharmaverseadamjnj

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

Generate ADaM datasets aligned with Johnson & Johnson's Clinical and Statistical Programming standards.


## Features

- Generates ADaM datasets that comply with J&J standards
- Built on top of the `pharmaverseadam` package
- Implements data conversion from pharmaverse format to J&J standards format
- Provides reproducible and consistent test data

## Implemented Datasets

Curently supports the following ADaM domains:
- ADSL (Subject Level Analysis Dataset)
- ADAE (Adverse Events Analysis Dataset)
- ADCM (Concomitant Medications Analysis Dataset)
- ADEG (ECG Analysis Dataset)
- ADEX (Exposure Analysis Dataset)
- ADESUM (Exposure Analysis Summary Dataset) 
- ADLB (Laboratory Test Results Analysis Dataset)
- ADVS (Vital Signs Analysis Dataset)
- ADTTESAF (Time-to-Event Safety Analysis Dataset)
- ADAEFMQ (Adverse Events Analysis Dataset FDA Medical Query)

## Installation

You can install the development version from GitHub:

```r
install.packages("pharmaverseadamjnj") #CRAN Release
#install.packages("devtools") #or for the rolling dev release
#devtools::install_github("johnsonandjohnson/pharmaverseadamjnj", ref = "dev")
```

## Usage

```r
library(pharmaverseadamjnj)

# Access the ADaM datasets directly
# These datasets are loaded lazily when the package is loaded
head(adsl)
head(adae)
head(adlb)
# ... and similarly for other domains
```

## Data Sources
Test datasets have been sourced from the [pharmaverseadam](https://github.com/pharmaverse/pharmaverseadam) package, which utilized the data from the [pharmaversesdtm](https://github.com/pharmaverse/pharmaversesdtm) package and the [CDISC pilot project](https://github.com/cdisc-org/sdtm-adam-pilot-project).
