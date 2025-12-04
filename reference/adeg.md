# adeg

adeg modified from pharmaverseadam

## Usage

``` r
adeg
```

## Format

A data frame with 13536 rows and 71 variables:

- DOMAIN:

  Domain Abbreviation

- USUBJID:

  Unique Subject Identifier

- EGSEQ:

  Sequence Number

- EGTESTCD:

  ECG Test or Examination Short Name

- EGTEST:

  ECG Test or Examination Name

- EGORRES:

  Result or Finding in Original Units

- EGORRESU:

  Original Units

- EGSTRESC:

  Character Result/Finding in Std Format

- EGSTRESN:

  Numeric Result/Finding in Standard Units

- EGSTRESU:

  Standard Units

- EGSTAT:

  Completion Status

- EGLOC:

  Lead Location Used for Measurement

- EGBLFL:

  Baseline Flag

- VISITNUM:

  Visit Number

- VISIT:

  Visit Name

- VISITDY:

  Planned Study Day of Visit

- EGDTC:

  Date/Time of ECG

- EGDY:

  Study Day of ECG

- EGTPT:

  Planned Time Point Name

- EGTPTNUM:

  Planned Time Point Number

- EGELTM:

  Planned Elapsed Time from Time Point Ref

- EGTPTREF:

  Time Point Reference

- ADTM:

  Analysis Datetime

- ATMF:

  Analysis Time Imputation Flag

- ADY:

  Analysis Relative Day

- PARAMCD:

  Parameter Code

- AVAL:

  Analysis Value

- AVALC:

  Analysis Value (C)

- ADT:

  Analysis Date

- ATPTN:

  Analysis Timepoint (N)

- ATPT:

  Analysis Timepoint

- AVISIT:

  Analysis Visit

- AVISITN:

  Analysis Visit (N)

- DTYPE:

  Derivation Type

- ONTRTFL:

  On Treatment Record Flag

- ANRLO:

  Analysis Normal Range Lower Limit

- ANRHI:

  Analysis Normal Range Upper Limit

- ANRIND:

  Analysis Reference Range Indicator

- BASETYPE:

  Baseline Type

- ABLFL:

  Baseline Record Flag

- BASEC:

  Baseline Value (C)

- CHG:

  Change from Baseline

- PCHG:

  Percent Change from Baseline

- ANL01FL:

  Analysis Flag 01-Analysis Value

- TRTP:

  Planned Treatment

- TRTA:

  Actual Treatment

- ASEQ:

  Analysis Sequence Number

- AVALCAT1:

  Analysis Value Category 1

- AVALCA1N:

  Analysis Value Category 1 (N)

- CHGCAT1:

  Change from Baseline Category 1

- CHGCAT1N:

  Change from Baseline Category 1 (N)

- PARAM:

  Parameter

- PARAMN:

  Parameter (N)

- TRTEMFL:

  Treatment Emergent Analysis Flag

- ANL02FL:

  Analysis Flag 02-By Visit Value

- ANL03FL:

  Analysis Flag 03-Maximum Value

- APOBLFL:

  Post-Baseline Record Flag

- CRIT1:

  Analysis Criterion 1

- CRIT1FL:

  Criterion 1 Evaluation Result Flag

- CRIT2:

  Analysis Criterion 2

- CRIT2FL:

  Criterion 2 Evaluation Result Flag

- EGCLSIG:

  Clinically Significant

- BASE:

  Baseline Value

- BNRIND:

  Baseline Reference Range Indicator

- BASECAT1:

  Baseline Category 1

- TRT01A:

  Actual Treatment for Period 01

- SAFFL:

  Safety Population Flag

- STUDYID:

  Study Identifier

- AGE:

  Age

- SEX:

  Sex

- RACE_DECODE:

  Race

## Source

data from pharmaverseadam.

## See also

[`adae`](adae.md) [`adaeocmq`](adaeocmq.md) [`adagocmq`](adagocmq.md)
[`adcm`](adcm.md) `adeg` [`adex`](adex.md) [`adexsum`](adexsum.md)
[`adlb`](adlb.md) [`adpc`](adpc.md) [`adsl`](adsl.md)
[`adttesaf`](adttesaf.md) [`advs`](advs.md)\# nolint

## Examples

``` r
head(data("adeg"))
#> [1] "adeg"
```
