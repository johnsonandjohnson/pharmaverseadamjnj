# pharmaverseadamjnj 0.0.6

### Added
- Added ADISHUM Dataset
- Added DDPCDTHC and NCTXSDT variables to ADSL
- Added AESEVN, EXSTDT, DOSEDT, ADATRES, and NABSTAT variables to ADAE
- Added EXSTDT, ADATRES, and NABSTAT variables to ADEX
- Added ADATRES and NABSTAT variables to ADEXSUM
- Added vars to ADSL: LASTCTDT, UNBLNDDT, IMFL, SAFEXRS, FASEXRS, PPREXRS, PKEXRES, and IMEXRES
- Added restord labels to xpt
- Change ONTRFl for ADVS to first dose within 30 days
- Added STRAT1D, STRAT2D to ADSL
- Added ABODSYS1, ABODSYS2, ADECOD1, and ADECOD2 to ADEX
- Added CRITx (1-8) to ADVS
- Added VSCLSIG and RACE to ADVS
- Added AOCTIFL, AESCAT, SMQzzNAM, and CQzzNAM, AERELSx, AEACNSx, RELGR1, AEDRGSx, DOSSxDY, DOSSxU, DOSSxON to ADAE
- Added derived variables `AESHOSPP` and `AESHOSPR`
- Added ADVS: RESP TEMP and PULSE PARAMCD
- Added variables `COHORT`, `GROUP`, and `EOTDT` to ADSL.
- Added variables `CMSTRTPT` and `CMSTTPT` to ADCM.
- Added AOCT01FL, AOCT02FL, AOCT03FL, AOCS01FL, AOCS02FL, AOCS03FL to ADAE
- Added ADAECOMP and ADSLCOMP datasets
- Added ADDISP dataset
- Added RACE variable from ADSL to ADEG dataset
- Added PKRFL variable to ADPC dataset

### Changed
- Modified ADSL death-related variable derivations (DTHFL, DTHDT, DTHCAUS, DDPCDTHC)
- Modified ADSL to add DCTDT label
- Modified ADEXSUM to add DURFUPD parameter
- Removed *_DECODE vars in ADSL, ADLB, ADVS, ADAE, ADCM, ADEG, ADEX
- Releveled TRT01P and TRT01A (now Xanomeline Low Dose -> Xanomeline High Dose -> Placebo)
- Modified "General disorders and administration site conditions" to "Gastrointestinal disorders" in ADEX
- Modified ADVS on ORTHYP, ORTHYPS and  ORTHYPD parameters for TSFVIT04
- Modified RESCRNFL and assigned 01-701-1240 to Subject refused to sign informed consent
- Modified ADVS on RESP parameter for gsfvit02a and gsfvit02b
- Modified ADVS on week 20 and week 24 AVISIT derivation
- Modified ADVS on AVISITN variable
- Modified ADVS on CRIT7/CRIT7FL CRIT8/CRIT8FL values of DIABP parameter
- Modified ADVS AVISIT, AVISITN and added ADY derivation
- Modified REGION1 to "Northern America"
- Added DCSREAS and DCTREAS levels
- Modified ADAE to update derivation rules for AEREL, AOCCFL, AOCCPFL, AOCCSFL, TRDISCFL, AESCAT and SMQxxNAM
- Modified ADAEOCMQ to update derivation rule for OCMQNAM and AEDECOD
- Modified ADAGOCMQ for AEDECOD merge logic, ACAT1, and used conversion factor GLUC, HBA1C for standard units to conventional units.
- Modified ADDISP to update labels for SxSDY and SxEDY variables.
- Fix adeg CRIT1FL and CRIT2FL factor order
- Modified ADEG for AVAL variableto include "ABNORMAL, CS", "ABNORMAL, NCS", "NORMAL"
- Modifed ADEG for PARAMCD convered 'msec' to 'ms'
- Modified ADPC for CRIT1 variable
- Modified ADEX to add ATPT, ATPTN, ADURC, ARSDSDO, ATDPRP, and  ATDPRPU variables
- Modified ADEX for APERIOD and updated ATPT values
- Modified ADEX to add AVAMT, AVAMTU, ADOSFRQP, ADOSDLY, ARSDOSD, ACDOSE, ACDOSU, ASTDY,
  ALL ANL flags ex. ANL01FL, ANL02FL, ANL03FL, ANL04FL, INTCAT, SKPCAT, DLYCAT variables.
- Updated PARAMN values in ADVS domain
- Added CRIT6, CRIT6FL, CRIT7 & CRIT7FL on ADILLI for [AST or ALT] PARCAT
- Added CRIT2, CRIT2FL, CRIT3, CRIT3FL on ADILLI for TBILI & ALP PARCAT
- Modified ADEXSUM to add TRTDURD parameter
- Modified ECOCCUR, ADOSE, ACDOSE, ASCHDOSE, AVAMT, ADOSU, ACDOSU. ASCHDOSU and AVAMTU variables
- Modified PARAM unit to 'mg' in ADEXSUM dataset
- Modified PARAM and PARAMCD to add CUMDOSS1, CUMDOSS2


## [0.0.5] - 2026-05-25
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



