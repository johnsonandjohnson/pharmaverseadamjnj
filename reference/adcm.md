# adcm

adcm modified from pharmaverseadam

## Usage

``` r
adcm
```

## Format

A data frame with 7276 rows and 63 variables:

- DOMAIN:

  Domain Abbreviation

- USUBJID:

  Unique Subject Identifier

- CMSEQ:

  Sequence Number

- CMSPID:

  Sponsor-Defined Identifier

- CMTRT:

  Reported Name of Drug, Med, or Therapy

- CMDECOD:

  Standardized Medication Name

- CMINDC:

  Indication

- CMCLAS:

  Medication Class

- CMDOSE:

  Dose per Administration

- CMDOSU:

  Dose Units

- CMDOSFRQ:

  Dosing Frequency per Interval

- CMROUTE:

  Route of Administration

- VISITNUM:

  Visit Number

- VISIT:

  Visit Name

- VISITDY:

  Planned Study Day of Visit

- CMDTC:

  Date/Time of Collection

- CMSTDTC:

  Start Date/Time of Medication

- CMENDTC:

  End Date/Time of Medication

- CMSTDY:

  Study Day of Start of Medication

- CMENDY:

  Study Day of End of Medication

- CMENRTPT:

  End Relative to Reference Time Point

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

- ONTRTFL:

  On Treatment Record Flag

- PREFL:

  Pre-treatment Flag

- FUPFL:

  Follow-up Flag

- ANL01FL:

  Analysis Flag 01

- AOCCPFL:

  1st Occurrence of Preferred Term Flag

- APHASE:

  Phase

- APHASEN:

  Description of Phase N

- TRTP:

  Planned Treatment

- TRTA:

  Actual Treatment

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

[`adae`](adae.md) [`adaeocmq`](adaeocmq.md) [`adagocmq`](adagocmq.md)
`adcm` [`adeg`](adeg.md) [`adex`](adex.md) [`adexsum`](adexsum.md)
[`adlb`](adlb.md) [`adpc`](adpc.md) [`adsl`](adsl.md)
[`adttesaf`](adttesaf.md) [`advs`](advs.md)\# nolint

## Examples

``` r
head(data("adcm"))
#> [1] "adcm"
```
