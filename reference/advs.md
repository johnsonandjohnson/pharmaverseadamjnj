# advs

advs modified from pharmaverseadam

## Usage

``` r
advs
```

## Format

A data frame with 40702 rows and 78 variables:

- DOMAIN:

  Domain Abbreviation

- USUBJID:

  Unique Subject Identifier

- VSSEQ:

  Sequence Number

- VSTESTCD:

  Vital Signs Test Short Name

- VSTEST:

  Vital Signs Test Name

- VSPOS:

  Vital Signs Position of Subject

- VSORRES:

  Result or Finding in Original Units

- VSORRESU:

  Original Units

- VSSTRESC:

  Character Result/Finding in Std Format

- VSSTRESN:

  Numeric Result/Finding in Standard Units

- VSSTRESU:

  Standard Units

- VSSTAT:

  Completion Status

- VSLOC:

  Location of Vital Signs Measurement

- VSBLFL:

  Baseline Flag

- VISITNUM:

  Visit Number

- VISIT:

  Visit Name

- VISITDY:

  Planned Study Day of Visit

- VSDTC:

  Date/Time of Measurements

- VSDY:

  Study Day of Vital Signs

- VSTPT:

  Planned Time Point Name

- VSTPTNUM:

  Planned Time Point Number

- VSELTM:

  Planned Elapsed Time from Time Point Ref

- VSTPTREF:

  Time Point Reference

- ADT:

  Analysis Date

- ADY:

  Analysis Relative Day

- PARAMCD:

  Parameter Code

- AVAL:

  Analysis Value

- ATPTN:

  Analysis Timepoint (N)

- ATPT:

  Analysis Time Point

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

- A1LO:

  Analysis Range 1 Lower Limit

- A1HI:

  Analysis Range 1 Upper Limit

- ANRIND:

  Analysis Reference Range Indicator

- BASETYPE:

  Baseline Type

- ABLFL:

  Baseline Record Flag

- CHG:

  Change from Baseline

- PCHG:

  Percent Change from Baseline

- ANL01FL:

  Analysis Flag 01

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

- PARAM:

  Parameter

- PARAMN:

  Parameter (N)

- AVALC:

  Analysis Value (C)

- ANL02FL:

  Analysis Flag 02-By Visit Value

- APOBLFL:

  Post-Baseline Record Flag

- BASE:

  Baseline Value

- BNRIND:

  Baseline Reference Range Indicator

- ADTM:

  Analysis Date/Time

- CRIT1:

  Analysis Criterion 1

- CRIT1FL:

  Criterion 1 Evaluation Result Flag

- CRIT2:

  Analysis Criterion 2

- CRIT2FL:

  Criterion 2 Evaluation Result Flag

- CRIT3:

  Analysis Criterion 3

- CRIT3FL:

  Criterion 3 Evaluation Result Flag

- ATOXDSCL:

  Analysis Toxicity Description Low

- ATOXDSCH:

  Analysis Toxicity Description High

- ATOXGRL:

  Analysis Toxicity Grade Low

- ATOXGRH:

  Analysis Toxicity Grade High

- ATOXGR:

  Analysis Toxicity Grade

- ANL06FL:

  Analysis Flag 06-Minimum Value

- ANL05FL:

  Analysis Flag 05-Worst Tox Grade High

- ANL04FL:

  Analysis Flag 04-Worst Value

- ANL03FL:

  Analysis Flag 03-Maximum Value

- TRTEMFL:

  Treatment Emergent Analysis Flag

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
[`adcm`](adcm.md) [`adeg`](adeg.md) [`adex`](adex.md)
[`adexsum`](adexsum.md) [`adlb`](adlb.md) [`adpc`](adpc.md)
[`adsl`](adsl.md) [`adttesaf`](adttesaf.md) `advs`\# nolint

## Examples

``` r
head(data("advs"))
#> [1] "advs"
```
