# adaeocmq

adae modified from pharmaverseadam to include Office of New Drugs (OND)
Custom Medical Queries

## Usage

``` r
adaeocmq
```

## Format

A data frame with 1972 rows and 86 variables:

- DOMAIN:

  Domain Abbreviation

- USUBJID:

  Unique Subject Identifier

- AESEQ:

  Sequence Number

- AESPID:

  Sponsor-Defined Identifier

- AETERM:

  Reported Term for the Adverse Event

- AELLT:

  Lowest Level Term

- AELLTCD:

  Lowest Level Term Code

- AEDECOD:

  AEDECOD

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

- AEBODSYS:

  Body System or Organ Class

- AEBDSYCD:

  Body System or Organ Class Code

- AESOC:

  Primary System Organ Class

- AESOCCD:

  Primary System Organ Class Code

- AESEV:

  Severity/Intensity

- AESER:

  Serious Event

- AEACN:

  Action Taken with Study Treatment

- AEREL:

  Causality

- AEOUT:

  Outcome of Adverse Event

- AESCAN:

  Involves Cancer

- AESCONG:

  Congenital Anomaly or Birth Defect

- AESDISAB:

  Persist or Signif Disability/Incapacity

- AESDTH:

  Results in Death

- AESHOSP:

  Requires or Prolongs Hospitalization

- AESLIFE:

  Is Life Threatening

- AESOD:

  Occurred with Overdose

- AEDTC:

  Date/Time of Collection

- AESTDTC:

  Start Date/Time of Adverse Event

- AEENDTC:

  End Date/Time of Adverse Event

- AESTDY:

  Study Day of Start of Adverse Event

- AEENDY:

  Study Day of End of Adverse Event

- ASTDTM:

  Analysis Start Date/Time

- ASTDTF:

  Analysis Start Date Imputation Flag

- ASTTMF:

  Analysis Start Time Imputation Flag

- AENDTM:

  Analysis End Date/Time

- AENDTF:

  Analysis End Date Imputation Flag

- AENTMF:

  Analysis End Time Imputation Flag

- ASTDT:

  Analysis Start Date

- AENDT:

  Analysis End Date

- ASTDY:

  Analysis Start Relative Day

- AENDY:

  Analysis End Relative Day

- ADURN:

  Analysis Duration (N)

- ADURU:

  Analysis Duration Units

- LDOSEDTM:

  End Date/Time of Last Dose

- ASEV:

  Analysis Severity/Intensity

- AREL:

  Analysis Causality

- TRTEMFL:

  Treatment Emergent Analysis Flag

- ASEVN:

  Analysis Severity/Intensity (N)

- AOCCIFL:

  1st Max Sev./Int. Occurrence Flag

- AETOXGR:

  Standard Toxicity Grade

- AETOXGRN:

  Standard Toxicity Grade (N)

- AEACN_DECODE:

  Action Taken with Study Treatment

- DOSEDY:

  Day of Study Drug

- DOSEU:

  Treatment Dose Units

- DOSEON:

  Treatment Dose at Record Start

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

- OCMQNAM:

  Custom Medical Query Name

- OCMQSOC:

  Custom Medical Query System Organ Class

- OCMQCLSS:

  Custom Medical Query Scope

- GENSPMFL:

  Gender Specific OCMQ Male Flag

- GENSPFFL:

  Gender Specific OCMQ Female Flag

## Source

data from adae from pharmaverseadam and, FDA_OCMQ_Consolidated_List.rds
and FDA_OCMQ_References.rds

## See also

[`adae`](adae.md) `adaeocmq` [`adagocmq`](adagocmq.md) [`adcm`](adcm.md)
[`adeg`](adeg.md) [`adex`](adex.md) [`adexsum`](adexsum.md)
[`adlb`](adlb.md) [`adpc`](adpc.md) [`adsl`](adsl.md)
[`adttesaf`](adttesaf.md) [`advs`](advs.md)\# nolint

## Examples

``` r
head(data("adaeocmq"))
#> [1] "adaeocmq"
```
