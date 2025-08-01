#' @title adcm
#'
#' @description adcm modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 7276 rows and 62 variables:
#' \describe{
#'  \item{DOMAIN}{Domain Abbreviation}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{CMSEQ}{Sequence Number}
#'  \item{CMSPID}{Sponsor-Defined Identifier}
#'  \item{CMTRT}{Reported Name of Drug, Med, or Therapy}
#'  \item{CMDECOD}{Standardized Medication Name}
#'  \item{CMINDC}{Indication}
#'  \item{CMCLAS}{Medication Class}
#'  \item{CMDOSE}{Dose per Administration}
#'  \item{CMDOSU}{Dose Units}
#'  \item{CMDOSFRQ}{Dosing Frequency per Interval}
#'  \item{CMROUTE}{Route of Administration}
#'  \item{VISITNUM}{Visit Number}
#'  \item{VISIT}{Visit Name}
#'  \item{VISITDY}{Planned Study Day of Visit}
#'  \item{CMDTC}{Date/Time of Collection}
#'  \item{CMSTDTC}{Start Date/Time of Medication}
#'  \item{CMENDTC}{End Date/Time of Medication}
#'  \item{CMSTDY}{Study Day of Start of Medication}
#'  \item{CMENDY}{Study Day of End of Medication}
#'  \item{CMENRTPT}{End Relative to Reference Time Point}
#'  \item{ASTDTM}{Analysis Start Date/Time}
#'  \item{ASTDTF}{Analysis Start Date Imputation Flag}
#'  \item{ASTTMF}{Analysis Start Time Imputation Flag}
#'  \item{AENDTM}{Analysis End Date/Time}
#'  \item{AENDTF}{Analysis End Date Imputation Flag}
#'  \item{AENTMF}{Analysis End Time Imputation Flag}
#'  \item{ASTDT}{Analysis Start Date}
#'  \item{AENDT}{Analysis End Date}
#'  \item{ASTDY}{Analysis Start Relative Day}
#'  \item{AENDY}{Analysis End Relative Day}
#'  \item{ADURN}{Analysis Duration (N)}
#'  \item{ADURU}{Analysis Duration Units}
#'  \item{ONTRTFL}{On Treatment Record Flag}
#'  \item{PREFL}{Pre-treatment Flag}
#'  \item{FUPFL}{Follow-up Flag}
#'  \item{ANL01FL}{Analysis Flag 01}
#'  \item{AOCCPFL}{1st Occurrence of Preferred Term Flag}
#'  \item{APHASE}{Phase}
#'  \item{APHASEN}{Description of Phase N}
#'  \item{TRTP}{Planned Treatment}
#'  \item{TRTA}{Actual Treatment}
#'  \item{CMLVL1}{Preferred ATC Text for ATC Level 1}
#'  \item{CMLVL2}{Preferred ATC Text for ATC Level 2}
#'  \item{CMLVL3}{Preferred ATC Text for ATC Level 3}
#'  \item{CMLVL4}{Preferred ATC Text for ATC Level 4}
#'  \item{CMBASPRF}{Base Preferred Term}
#'  \item{CMPRESP}{CM Pre-specified}
#'  \item{CMOCCUR}{CM Occurrence}
#'  \item{CMINDCSP}{Indication Specification}
#'  \item{CMDOSTXT}{Dose Description}
#'  \item{CMENRF}{End Relative to Reference Period}
#'  \item{CQ01NAM}{Customized Query 01 Name}
#'  \item{CQ02NAM}{Customized Query 02 Name}
#'  \item{CQ03NAM}{Customized Query 03 Name}
#'  \item{CQ04NAM}{Customized Query 04 Name}
#'  \item{CQ05NAM}{Customized Query 05 Name}
#'  \item{CQ06NAM}{Customized Query 06 Name}
#'  \item{CQ07NAM}{Customized Query 07 Name}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{TRTSDT}{Date of First Exposure to Treatment}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaefmq}} \code{\link{adcm}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets adcm
#' @name adcm
#' @examples
#'  head(data("adcm"))
"adcm"
