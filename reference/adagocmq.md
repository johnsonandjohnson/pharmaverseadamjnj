# adagocmq

adagocmq modified from pharmaverseadam

## Usage

``` r
adagocmq
```

## Format

A data frame with 1301 rows and 27 variables:

- STUDYID:

  Study Identifier

- USUBJID:

  Unique Subject Identifier

- ASTDT:

  Analysis Start Date

- ASTDY:

  Analysis Start Relative Day

- TRT01P:

  Planned Treatment for Period 01

- TRT01PN:

  Planned Treatment for Period 01 (N)

- TRT01A:

  Actual Treatment for Period 01

- TRT01AN:

  Actual Treatment for Period 01 (N)

- ACAT1:

  Analysis Category 1

- ACAT1N:

  Analysis Category 1 (N)

- ATERM:

  Analysis Term (N)

- ATERMN:

  Analysis Term

- HYPSCAT:

  Hypersensitivity Category

- SRCVALUE:

  Source Value

- SRCVAR:

  Source Variable

- SRCDOM:

  Source Data

- SRCSEQ:

  Source Sequence Number

- ANL01FL:

  Analysis Flag 01

- AGE:

  Age

- AGEU:

  Age Units

- SEX:

  Sex

- RACE:

  Race

- COUNTRY:

  Country

- RANDFL:

  Randomized Flag

- SAFFL:

  Safety Population Flag

- SITEID:

  Study Site Identifier

- SUBJID:

  Subject Identifier for the Study

## Source

data from pharmaverseadam.

## See also

[`adae`](adae.md) [`adaeocmq`](adaeocmq.md) `adagocmq` [`adcm`](adcm.md)
[`adeg`](adeg.md) [`adex`](adex.md) [`adexsum`](adexsum.md)
[`adlb`](adlb.md) [`adpc`](adpc.md) [`adsl`](adsl.md)
[`adttesaf`](adttesaf.md) [`advs`](advs.md)\# nolint

## Examples

``` r
head(data("adagocmq"))
#> [1] "adagocmq"
```
