# adcm

adcm modified from pharmaverseadam

## Usage

``` r
adcm
```

## Format

A data frame with 7276 rows and 63 variables:

- USUBJID:

  Unique Subject Identifier

- DOMAIN:

  Domain Abbreviation

- TRTP:

  Planned Treatment

- TRTA:

  Actual Treatment

- APHASE:

  Phase

- APHASEN:

  Description of Phase N

- CMSEQ:

  Sequence Number

- CMDECOD:

  Standardized Medication Name

- CMTRT:

  Reported Name of Drug, Med, or Therapy

- CMCLAS:

  Medication Class

- CMSTDTC:

  Start Date/Time of Medication

- ASTDT:

  Analysis Start Date

- ASTDTM:

  Analysis Start Date/Time

- ASTDTF:

  Analysis Start Date Imputation Flag

- ASTTMF:

  Analysis Start Time Imputation Flag

- CMENDTC:

  End Date/Time of Medication

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

- CMSTDY:

  Study Day of Start of Medication

- AENDY:

  Analysis End Relative Day

- CMENDY:

  Study Day of End of Medication

- ADURN:

  Analysis Duration (N)

- ADURU:

  Analysis Duration Units

- ANL01FL:

  Analysis Flag 01

- ONTRTFL:

  On Treatment Record Flag

- PREFL:

  Pre-treatment Flag

- FUPFL:

  Follow-up Flag

- AOCCPFL:

  1st Occurrence of Preferred Term Flag

- CMINDC:

  Indication

- CMDOSE:

  Dose per Administration

- CMDOSU:

  Dose Units

- CMDOSFRQ:

  Dosing Frequency per Interval

- CMROUTE:

  Route of Administration

- CMSPID:

  Sponsor-Defined Identifier

- CMENRTPT:

  End Relative to Reference Time Point

- VISITNUM:

  Visit Number

- VISIT:

  Visit Name

- VISITDY:

  Planned Study Day of Visit

- CMDTC:

  Date/Time of Collection

- CMLVL1:

  Preferred ATC Text for ATC Level 1

- CMLVL2:

  Preferred ATC Text for ATC Level 2

- CMLVL3:

  Preferred ATC Text for ATC Level 3

- CMLVL4:

  Preferred ATC Text for ATC Level 4

- CMBASPRF:

  Base Preferred Term

- CMPRESP:

  CM Pre-specified

- CMOCCUR:

  CM Occurrence

- CMINDCSP:

  Indication Specification

- CMDOSTXT:

  Dose Description

- CMENRF:

  End Relative to Reference Period

- CQ01NAM:

  Customized Query 01 Name

- CQ02NAM:

  Customized Query 02 Name

- CQ03NAM:

  Customized Query 03 Name

- CQ04NAM:

  Customized Query 04 Name

- CQ05NAM:

  Customized Query 05 Name

- CQ06NAM:

  Customized Query 06 Name

- CQ07NAM:

  Customized Query 07 Name

- TRT01A:

  Actual Treatment for Period 01

- SAFFL:

  Safety Population Flag

- TRTSDT:

  Date of First Exposure to Treatment

- STUDYID:

  Study Identifier

## Source

data from pharmaverseadam.

## See also

[`adae`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adae.md)
[`adaeocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adaeocmq.md)
[`adagocmq`](https://johnsonandjohnson.github.io/pharmaverseadamjnj/reference/adagocmq.md)
`adcm`
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
head(data("adcm"))
#> [1] "adcm"
```
