# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a
Changelog](https://keepachangelog.com/en/1.0.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## \[0.0.2\] - 2025-12-03

### Added

- Added ADAGOCMQ & ADPC dataset
- Added `LBSPEC` and `LBFAST` to ADLB dataset

### Changed

- Update ADSL to have all factor levels of RACE and ETHNIC
- Update adsl DTH60TFL to DTHB60FL
- Updated adae update to include start/end time.
- Updated source file OCMQs to reflect FDA changes to FMQ
- Renamed adaefmq to adaeocmq
- Fixed ADCM to properly handle uncoded terms
- Fixed ADLB to use “BEFORE TREATMENT” for ATPT instead of “00:00”
- ADSL: Set RANDFL to “Y” for most subjects with SAFFL=“Y”, but also for
  a few with SAFFL=“N” and non-NA TRT01P

## \[0.0.1\] - 2024-07-02

### Added

- Initial CRAN release
- Provides ADaM datasets that comply with J&J Innovative Medicine
  standards
- Built on top of the ‘pharmaverseadam’ package
- Implements data conversion from pharmaverse format to J&J standards
  format
- Provides reproducible and consistent test data

### Implemented ADaM Domains

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
