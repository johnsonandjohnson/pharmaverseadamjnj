# adeg

adeg modified from pharmaverseadam

## Usage

``` r
adeg
```

## Format

A data frame with 13536 rows and 71 variables:

- USUBJID:

  Unique Subject Identifier

- DOMAIN:

  Domain Abbreviation

- ASEQ:

  Analysis Sequence Number

- TRTP:

  Planned Treatment

- TRTA:

  Actual Treatment

- ADT:

  Analysis Date

- ADTM:

  Analysis Datetime

- ADY:

  Analysis Relative Day

- ATMF:

  Analysis Time Imputation Flag

- AVISIT:

  Analysis Visit

- AVISITN:

  Analysis Visit (N)

- ATPT:

  Analysis Timepoint

- ATPTN:

  Analysis Timepoint (N)

- PARAM:

  Parameter

- PARAMCD:

  Parameter Code

- PARAMN:

  Parameter (N)

- AVAL:

  Analysis Value

- AVALC:

  Analysis Value (C)

- AVALCAT1:

  Analysis Value Category 1

- AVALCA1N:

  Analysis Value Category 1 (N)

- BASEC:

  Baseline Value (C)

- BASETYPE:

  Baseline Type

- CHG:

  Change from Baseline

- CHGCAT1:

  Change from Baseline Category 1

- CHGCAT1N:

  Change from Baseline Category 1 (N)

- PCHG:

  Percent Change from Baseline

- DTYPE:

  Derivation Type

- ANRIND:

  Analysis Reference Range Indicator

- ANRLO:

  Analysis Normal Range Lower Limit

- ANRHI:

  Analysis Normal Range Upper Limit

- ABLFL:

  Baseline Record Flag

- ANL01FL:

  Analysis Flag 01-Analysis Value

- ONTRTFL:

  On Treatment Record Flag

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

[`adae`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adae.md)
[`adaeocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adaeocmq.md)
[`adagocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adagocmq.md)
[`adcm`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adcm.md)
[`addili`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/addili.md)
`adeg`
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
head(data("adeg"))
#> [1] "adeg"
```
