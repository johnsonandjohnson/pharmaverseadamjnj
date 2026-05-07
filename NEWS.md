# pharmaverseadamjnj 0.0.4

### Added
- Added vars to ADSL: LASTCTDT, UNBLNDDT, IMFL, SAFEXRS, FASEXRS, PPREXRS, PKEXRES, and IMEXRES
- Added restord labels to xpt
- Change ONTRFl for ADVS to first dose within 30 days
- Added STRAT1D, STRAT2D to ADSL
- Added ABODSYS1, ABODSYS2, ADECOD1, and ADECOD2 to ADEX
- Added CRITx (1-8) to ADVS
- Added VSCLSIG and RACE to ADVS
- Added AOCTIFL, AESCAT, SMQzzNAM, and CQzzNAM to ADAE
- Added derived variables `AESHOSPP` and `AESHOSPR`
- Added ADVS: RESP TEMP and PULSE PARAMCD
- Added variables `COHORT`, `GROUP`, and `EOTDT` to ADSL.
- Added variables `CMSTRTPT` and `CMSTTPT` to ADCM.


### Changed
- Removed *_DECODE vars in ADSL, ADLB, ADVS, ADAE, ADCM, ADEG, ADEX
- Releveled TRT01P and TRT01A (now Xanomeline Low Dose -> Xanomeline High Dose -> Placebo)
- Modified "General disorders and administration site conditions" to "Gastrointestinal disorders" in ADEX
- Modified ADVS on ORTHYP, ORTHYPS and  ORTHYPD parameters for TSFVIT04
- Modified RESCRNFL and assigned 01-701-1240 to Subject refused to sign informed consent

## [0.0.3] - 2026-01-27

### Added
- Add DCTADY to ADSL #32
- Add new addili dataset
- Add LBNAM variable to adlb 
- Add REGION1 = "North America"


## [0.0.2] - 2025-12-03

- Added ADAGOCMQ &  ADPC dataset
- Added `LBSPEC` and `LBFAST` to ADLB dataset

### Changed
- Update ADSL to have all factor levels of RACE and ETHNIC
- Update adsl DTH60TFL to DTHB60FL
- Updated adae update to include start/end time.
- Updated source file OCMQs to reflect FDA changes to FMQ
- Renamed adaefmq to adaeocmq
- Fixed ADCM to properly handle uncoded terms
- Fixed ADLB to use "BEFORE TREATMENT" for ATPT instead of "00:00"
- ADSL: Set RANDFL to "Y" for most subjects with SAFFL="Y", but also for a few with SAFFL="N" and non-NA TRT01P
- Updated {pharmaverseadam} released a new version 0.1.2 some variable labels have been added accordingly.

## [0.0.1] - 2024-07-02

### Added
- Initial CRAN release
- Provides ADaM datasets that comply with J&J Innovative Medicine standards
- Built on top of the 'pharmaverseadam' package
- Implements data conversion from pharmaverse format to J&J standards format
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



