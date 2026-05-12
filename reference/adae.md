# adae

adae modified from pharmaverseadam

## Usage

``` r
adae
```

## Format

A data frame with 1191 rows and 81 variables:

- USUBJID:

  Unique Subject Identifier

- DOMAIN:

  Domain Abbreviation

- AESEQ:

  Sequence Number

- AETERM:

  Reported Term for the Adverse Event

- AEDECOD:

  Dictionary-Derived Term

- AEBODSYS:

  Body System or Organ Class

- AEBDSYCD:

  Body System or Organ Class Code

- AELLT:

  Lowest Level Term

- AELLTCD:

  Lowest Level Term Code

- AEPTCD:

  Preferred Term Code

- AEHLT:

  High Level Term

- AEHLTCD:

  High Level Term Code

- AEHLGT:

  High Level Group Term

- AEHLGTCD:

  High Level Group Term Code

- AESOC:

  Primary System Organ Class

- AESOCCD:

  Primary System Organ Class Code

- AESTDTC:

  Start Date/Time of Adverse Event

- ASTDT:

  Analysis Start Date

- ASTDTM:

  Analysis Start Date/Time

- ASTDTF:

  Analysis Start Date Imputation Flag

- ASTTMF:

  Analysis Start Time Imputation Flag

- AEENDTC:

  End Date/Time of Adverse Event

- AENDT:

  Analysis End Date

- AENDTM:

  Analysis End Date/Time

- AENDTF:

  Analysis End Date Imputation Flag

- AENTMF:

  Analysis End Time Imputation Flag

- ASTDY:

  Analysis Start Relative Day

- AESTDY:

  Study Day of Start of Adverse Event

- AENDY:

  Analysis End Relative Day

- AEENDY:

  Study Day of End of Adverse Event

- ADURN:

  Analysis Duration (N)

- ADURU:

  Analysis Duration Units

- TRTEMFL:

  Treatment Emergent Analysis Flag

- AOCCIFL:

  1st Max Sev./Int. Occurrence Flag

- AESER:

  Serious Event

- AESDTH:

  Results in Death

- AESLIFE:

  Is Life Threatening

- AESHOSP:

  Requires or Prolongs Hospitalization

- AESDISAB:

  Persist or Signif Disability/Incapacity

- AESCONG:

  Congenital Anomaly or Birth Defect

- AESEV:

  Severity/Intensity

- ASEV:

  Analysis Severity/Intensity

- ASEVN:

  Analysis Severity/Intensity (N)

- AEREL:

  Causality

- AREL:

  Analysis Causality

- AEACN:

  Action Taken with Study Treatment

- AESPID:

  Sponsor-Defined Identifier

- AEOUT:

  Outcome of Adverse Event

- AESCAN:

  Involves Cancer

- AESOD:

  Occurred with Overdose

- AEDTC:

  Date/Time of Collection

- LDOSEDTM:

  End Date/Time of Last Dose

- DOSEON:

  Treatment Dose at Record Start

- DOSEU:

  Treatment Dose Units

- AETOXGR:

  Standard Toxicity Grade

- AETOXGRN:

  Standard Toxicity Grade (N)

- AEACN_DECODE:

  Action Taken with Study Treatment

- DOSEDY:

  Day of Study Drug

- AECONTRT:

  Concomitant or Additional Trtmnt Given

- CQ01NAM:

  Customized Query 01 Name

- CQ02NAM:

  Customized Query 02 Name

- CQ03NAM:

  Customized Query 03 Name

- AESMIE:

  Other Medically Important Serious Event

- AESER_DECODE:

  Serious Event

- AEREL_DECODE:

  Causality

- AEOUT_DECODE:

  Outcome of Adverse Event

- ACAT1:

  Analysis Category 1

- AOCCFL:

  1st Occurrence within Subject Flag

- AOCCPFL:

  1st Occurrence within Preferred Term Flag

- AOCCSFL:

  1st Occurrence of SOC Flag

- TRT01A:

  Actual Treatment for Period 01

- SAFFL:

  Safety Population Flag

- AGE:

  Age

- SEX:

  Sex

- RACE:

  Race

- RACE_DECODE:

  Race

- STUDYID:

  Study Identifier

- AGEGR1:

  Pooled Age Group 1

- TRTEDY:

  Treatment Relative End Day

- TRT01P:

  Planned Treatment for Period 01

- TRDISCFL:

  Treatment Discontinued Flag

## Source

data from pharmaverseadam.

## See also

`adae`
[`adaeocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adaeocmq.md)
[`adagocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adagocmq.md)
[`adcm`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adcm.md)
[`addili`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/addili.md)
[`adeg`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adeg.md)
[`adex`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adex.md)
[`adexsum`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adexsum.md)
[`adlb`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adlb.md)
[`adpc`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adpc.md)
[`adsl`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adsl.md)
[`adttesaf`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adttesaf.md)
[`advs`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/advs.md)\#
nolint

## Examples

``` r
head(data("adae"))
#> [1] "adae"
```
