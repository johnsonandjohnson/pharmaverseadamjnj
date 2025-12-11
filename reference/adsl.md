# adsl

adsl modified from pharmaverseadam

## Usage

``` r
adsl
```

## Format

A data frame with 306 rows and 108 variables:

- STUDYID:

  Study Identifier

- USUBJID:

  Unique Subject Identifier

- SUBJID:

  Subject Identifier for the Study

- RFSTDTC:

  Subject Reference Start Date/Time

- RFENDTC:

  Subject Reference End Date/Time

- RFXSTDTC:

  Date/Time of First Study Treatment

- RFXENDTC:

  Date/Time of Last Study Treatment

- RFICDTC:

  Date/Time of Informed Consent

- RFPENDTC:

  Date/Time of End of Participation

- DTHDTC:

  Date/Time of Death

- DTHFL:

  Subject Death Flag

- SITEID:

  Study Site Identifier

- AGE:

  Age

- AGEU:

  Age Units

- SEX:

  Sex

- RACE:

  Race

- ETHNIC:

  Ethnicity

- ARMCD:

  Planned Arm Code

- ARM:

  Description of Planned Arm

- ACTARMCD:

  Actual Arm Code

- ACTARM:

  Description of Actual Arm

- COUNTRY:

  Country

- DMDTC:

  Date/Time of Collection

- DMDY:

  Study Day of Collection

- TRT01P:

  Planned Treatment for Period 01

- TRT01A:

  Actual Treatment for Period 01

- TRTSDTM:

  Datetime of First Exposure to Treatment

- TRTSTMF:

  Time of First Exposure Imput. Flag

- TRTEDTM:

  Datetime of Last Exposure to Treatment

- TRTETMF:

  Time of Last Exposure Imput. Flag

- TRTSDT:

  Date of First Exposure to Treatment

- TRTEDT:

  Date of Last Exposure to Treatment

- TRTDURD:

  Total Treatment Duration (Days)

- SCRFDT:

  Screen Failure Date

- EOSDT:

  End of Study Date

- EOSSTT:

  End of Study Status

- FRVDT:

  Final Retrieval Visit Date

- RANDDT:

  Date of Randomization

- DTHDT:

  Date of Death

- DTHDTF:

  Date of Death Imputation Flag

- DTHADY:

  Relative Day of Death

- LDDTHELD:

  Elapsed Days from Last Dose to Death

- DTHCAUS:

  Cause of Death

- DTHDOM:

  Domain for Date of Death Collection

- DTHCGR1:

  Cause of Death Reason 1

- LSTALVDT:

  Date Last Known Alive

- SAFFL:

  Safety Population Flag

- RACEGR1:

  Pooled Race Group 1

- AGEGR1:

  Pooled Age Group 1

- REGION1:

  Geographic Region 1

- LDDTHGR1:

  Last Dose to Death - Days Elapsed Grp 1

- DTH30FL:

  Death Within 30 Days of Last Trt Flag

- DTHA30FL:

  Death After 30 Days from Last Trt Flag

- DTHB30FL:

  Death Within 30 Days of First Trt Flag

- TRT01PN:

  Planned Treatment for Period 01 (N)

- TRT01AN:

  Actual Treatment for Period 01 (N)

- AGEGR1N:

  Pooled Age Group 1 (N)

- SEX_DECODE:

  Sex

- WEIGHTBL:

  Weight (kg)

- WGTGR1N:

  Weight Group 1 (N)

- WGTGR1:

  Weight Group 1

- HEIGHTBL:

  Height (cm)

- BSABL:

  Body surface area (m2)

- BMIBL:

  Body mass index (kg/m2)

- BMIBLG1N:

  BMI at Baseline Group 1 (N)

- BMIBLG1:

  BMI at Baseline Group 1

- COUNTRY_DECODE:

  Country

- RACE_DECODE:

  Race

- RFICDT:

  Date of Informed Consent

- ETHNIC_DECODE:

  Ethnicity

- STRAT1R:

  Strat Factor 1 Value Used for Rand

- STRAT2R:

  Strat Factor 2 Value Used for Rand

- RANUM:

  Randomization Number

- RANDDTM:

  Datetime of Randomization

- EOTSTT:

  End of Treatment Status

- DCTREAS:

  Reason for Discontinuation of Treatment

- LTVISIT:

  Last Treatment Visit

- DCTREASP:

  Reason Specify for Discont of Treatment

- DCTDT:

  End of Study Date

- DCSREAS:

  Reason for Discontinuation from Study

- DCSREASP:

  Reason Spec for Discont from Study

- LSVISIT:

  Last Study Visit

- TRTEDY:

  Treatment Relative End Day

- SCRNFL:

  Screened Population Flag

- SCRFFL:

  Screen Failure Flag

- DCSCREEN:

  Reason for Discont During Screening

- ENRLFL:

  Enrolled Population Flag

- RANDFL:

  Randomized Flag

- ITTFL:

  Intent-To-Treat Population Flag

- FASFL:

  Full Analysis Set Population Flag

- PPROTFL:

  Per-Protocol Population Flag

- LSTSVDT:

  Last Subject Visit (SV) Date

- EOSDY:

  Study Day of Study Termination

- UNBLNDFL:

  Subject Blind Broken

- RESCRNFL:

  Re-screened Flag

- DTHTRTFL:

  Death on Treatment Flag

- DTHCAUSP:

  Cause Spec for Death

- DTHAFTFL:

  Death After 30 Days of Last Treatment

- DTHB60FL:

  Death Within 60 Days of First Treatment

- UNBLNDDY:

  Study Day of Unblinding

- UNBREAS:

  Reason For Unblinding

- LDOSE:

  Last Dose

- LDOSU:

  Last Dose Unit

- DTHTERM:

  Reported Cause of Death

- LDSTODTH:

  Days from Last Dose to Death

- DTHDY:

  Study Day of Death

- PKFL:

  Pharmacokinetic Population Flag

- DIABETFL:

  History of Diabetes

## Source

data from pharmaverseadam.

## See also

[`adae`](adae.md) [`adaeocmq`](adaeocmq.md) [`adagocmq`](adagocmq.md)
[`adcm`](adcm.md) [`adeg`](adeg.md) [`adex`](adex.md)
[`adexsum`](adexsum.md) [`adlb`](adlb.md) [`adpc`](adpc.md) `adsl`
[`adttesaf`](adttesaf.md) [`advs`](advs.md)\# nolint

## Examples

``` r
head(data("adsl"))
#> [1] "adsl"
```
