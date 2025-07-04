# pharmaverseadamjnj

<!-- start badges -->
[![Check 🛠](https://github.com/johnsonandjohnson/pharmaverseadamjnj/actions/workflows/check.yaml/badge.svg)](https://github.com/johnsonandjohnson/pharmaverseadamjnj/actions/workflows/check.yaml)
[![Docs 📚](https://github.com/johnsonandjohnson/pharmaverseadamjnj/actions/workflows/pkgdown.yaml/badge.svg)](https://johnsonandjohnson.github.io/pharmaverseadamjnj/)

![GitHub forks](https://img.shields.io/github/forks/johnsonandjohnson/pharmaverseadamjnj?style=social)
![GitHub repo stars](https://img.shields.io/github/stars/johnsonandjohnson/pharmaverseadamjnj?style=social)

![GitHub commit activity](https://img.shields.io/github/commit-activity/m/johnsonandjohnson/pharmaverseadamjnj)
![GitHub contributors](https://img.shields.io/github/contributors/johnsonandjohnson/pharmaverseadamjnj)
![GitHub last commit](https://img.shields.io/github/last-commit/johnsonandjohnson/pharmaverseadamjnj)
![GitHub pull requests](https://img.shields.io/github/issues-pr/johnsonandjohnson/pharmaverseadamjnj)
![GitHub repo size](https://img.shields.io/github/repo-size/johnsonandjohnson/pharmaverseadamjnj)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Current Version](https://img.shields.io/github/r-package/v/johnsonandjohnson/pharmaverseadamjnj/main?color=purple&label=package%20version)](https://github.com/johnsonandjohnson/pharmaverseadamjnj/tree/main)
[![Open Issues](https://img.shields.io/github/issues-raw/johnsonandjohnson/pharmaverseadamjnj?color=red&label=open%20issues)](https://github.com/johnsonandjohnson/pharmaverseadamjnj/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc)
<!-- [![Coverage](https://github.com/johnsonandjohnson/pharmaverseadamjnj/actions/workflows/coverage.yaml/badge.svg)](https://github.com/johnsonandjohnson/pharmaverseadamjnj/actions/workflows/coverage.yaml) -->
<!-- end badges -->

Generate ADaM datasets aligned with Johnson & Johnson's Clinical and Statistical Programming standards.


## Features

- Generates ADaM datasets that comply with J&J standards
- Built on top of the 'pharmaverseadam' package
- Implements data conversion from pharmaverse format to J&J standards format
- Provides reproducible and consistent test data

## Implemented Datasets

Currently supports the following ADaM domains:

* ADSL (Subject Level Analysis Dataset)
* ADAE (Adverse Events Analysis Dataset)
* ADCM (Concomitant Medications Analysis Dataset)
* ADEG (ECG Analysis Dataset)
* ADEX (Exposure Analysis Dataset)
* ADESUM (Exposure Analysis Summary Dataset) 
* ADLB (Laboratory Test Results Analysis Dataset)
* ADVS (Vital Signs Analysis Dataset)
* ADTTESAF (Time-to-Event Safety Analysis Dataset)
* ADAEFMQ (Adverse Events Analysis Dataset FDA Medical Query)

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("johnsonandjohnson/pharmaverseadamjnj")
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
