# adishum

adishum modified from pharmaverseadam

## Usage

``` r
adishum
```

## Format

A data frame with 10922 rows and 66 variables:

- USUBJID:

  Unique Subject Identifier

- VISITNUM:

  Visit Number

- PARAMCD:

  Parameter Code

- ISTESTCD:

  Immunogenicity Test/Exam Short Name

- PARCAT1:

  Parameter Category 1

- PARAM:

  Parameter Description

- ISSEQ:

  Sequence Number

- ISBDAGNT:

  Binding Agent

- ISSTRESC:

  Character Result/Finding in Std Format

- ISSTRESN:

  Numeric Results/Findings in Std. Units

- ISNAM:

  Vendor Name

- VISIT:

  Visit Name

- ISDTC:

  Date/Time of Collection

- ADTM:

  Analysis Datetime

- PARQUAL:

  Parameter Qualifier

- AVAL:

  Analysis Value

- AVALC:

  Analysis Value (C)

- AVALCAT1:

  Analysis Value Category 1

- IMEVFL:

  IS Evaluable Population Flag

- AVISITN:

  Analysis Visit (N)

- AVISIT:

  Analysis Visit

- VISITDY:

  Planned Study Day of Visit

- ADT:

  Analysis Date

- ADY:

  Analysis Relative Day

- ABLFL:

  Baseline Record Flag

- ANL01FL:

  Analysis Flag 01

- ANL02FL:

  Analysis Flag 02

- ANL03FL:

  Analysis Flag 03-Inf SR

- ANL04FL:

  Analysis Flag 04-Severe Inf SR

- ANL05FL:

  Analysis Flag 05-Serious Inf SR

- ANL06FL:

  Analysis Flag 06-Infus SR - DC

- ANL07FL:

  Analysis Flag 07-Inj SR

- ANL08FL:

  Analysis Flag 08-Severe Inj SR

- ANL09FL:

  Analysis Flag 09-Serious Inj SR

- ANL10FL:

  Analysis Flag 10-Inject SR - DC

- ANL11FL:

  Analysis Flag 11-grade \>=3 Inf SR

- ANL12FL:

  Analysis Flag 12-grade \>=3 Inj SR

- ATPT:

  Analysis Timepoint

- ATPTN:

  Analysis Timepoint (N)

- INFUSPN:

  Number of Placebo infusions

- INFUSPRN:

  Placebo infusions with reaction

- INFUSAN:

  Number of Active Drug infusions

- INFUSARN:

  Active Drug infusions with reaction

- INJECPN:

  Number of Placebo injections

- INJECPRN:

  Placebo injections with reaction

- INJECAN:

  Number of Active Drug injections

- INJECARN:

  Active Drug injections with reaction

- DTL:

  Drug Tolerance Limit

- ADATRES:

  Treatment-emergent ADA Subject Status

- NABSTAT:

  NAB Status

- ADATREPT:

  TE ADA Positive Subject Peak Titer

- STUDYID:

  Study Identifier

- SUBJID:

  Subject Identifier for the Study

- SITEID:

  Study Site Identifier

- AGE:

  Age

- SEX:

  Sex

- RACE:

  Race

- IMFL:

  Immunogenicity Population Flag

- SAFFL:

  Safety Population Flag

- TRTSDT:

  Date of First Exposure to Treatment

- TRTSDTM:

  Datetime of First Exposure to Treatment

- TRT01P:

  Planned Treatment for Period 01

- TRT01PN:

  Planned Treatment for Period 01 (N)

- TRT01A:

  Actual Treatment for Period 01

- TRT01AN:

  Actual Treatment for Period 01 (N)

- TRTA:

  Actual Treatment

## Source

data from pharmaverseadam.

## See also

[`adae`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adae.md)
[`adaecomp`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adaecomp.md)
[`adaeocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adaeocmq.md)
[`adagocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adagocmq.md)
[`adcm`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adcm.md)
[`addili`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/addili.md)
[`addisp`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/addisp.md)
[`adeg`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adeg.md)
[`adex`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adex.md)
[`adexsum`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adexsum.md)
`adishum`
[`adlb`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adlb.md)
[`adpc`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adpc.md)
[`adsl`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adsl.md)
[`adslcomp`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adslcomp.md)
[`adttesaf`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adttesaf.md)
[`advs`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/advs.md)\#
nolint

## Examples

``` r
 head(data("adishum"))
#> [1] "adishum"
```
