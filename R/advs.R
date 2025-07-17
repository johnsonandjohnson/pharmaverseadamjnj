#' @title advs
#'
#' @description advs modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 40702 rows and 78 variables:
#' \describe{
#'  \item{DOMAIN}{Domain Abbreviation}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{VSSEQ}{Sequence Number}
#'  \item{VSTESTCD}{Vital Signs Test Short Name}
#'  \item{VSTEST}{Vital Signs Test Name}
#'  \item{VSPOS}{Vital Signs Position of Subject}
#'  \item{VSORRES}{Result or Finding in Original Units}
#'  \item{VSORRESU}{Original Units}
#'  \item{VSSTRESC}{Character Result/Finding in Std Format}
#'  \item{VSSTRESN}{Numeric Result/Finding in Standard Units}
#'  \item{VSSTRESU}{Standard Units}
#'  \item{VSSTAT}{Completion Status}
#'  \item{VSLOC}{Location of Vital Signs Measurement}
#'  \item{VSBLFL}{Baseline Flag}
#'  \item{VISITNUM}{Visit Number}
#'  \item{VISIT}{Visit Name}
#'  \item{VISITDY}{Planned Study Day of Visit}
#'  \item{VSDTC}{Date/Time of Measurements}
#'  \item{VSDY}{Study Day of Vital Signs}
#'  \item{VSTPT}{Planned Time Point Name}
#'  \item{VSTPTNUM}{Planned Time Point Number}
#'  \item{VSELTM}{Planned Elapsed Time from Time Point Ref}
#'  \item{VSTPTREF}{Time Point Reference}
#'  \item{ADT}{Analysis Date}
#'  \item{ADY}{Analysis Relative Day}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{AVAL}{Analysis Value}
#'  \item{ATPTN}{Analysis Timepoint (N)}
#'  \item{ATPT}{Analysis Time Point}
#'  \item{AVISIT}{Analysis Visit}
#'  \item{AVISITN}{Analysis Visit (N)}
#'  \item{DTYPE}{Derivation Type}
#'  \item{ONTRTFL}{On Treatment Record Flag}
#'  \item{ANRLO}{Analysis Normal Range Lower Limit}
#'  \item{ANRHI}{Analysis Normal Range Upper Limit}
#'  \item{A1LO}{Analysis Range 1 Lower Limit}
#'  \item{A1HI}{Analysis Range 1 Upper Limit}
#'  \item{ANRIND}{Analysis Reference Range Indicator}
#'  \item{BASETYPE}{Baseline Type}
#'  \item{ABLFL}{Baseline Record Flag}
#'  \item{CHG}{Change from Baseline}
#'  \item{PCHG}{Percent Change from Baseline}
#'  \item{ANL01FL}{Analysis Flag 01}
#'  \item{TRTP}{Planned Treatment}
#'  \item{TRTA}{Actual Treatment}
#'  \item{ASEQ}{Analysis Sequence Number}
#'  \item{AVALCAT1}{Analysis Value Category 1}
#'  \item{AVALCA1N}{Analysis Value Category 1 (N)}
#'  \item{PARAM}{Parameter}
#'  \item{PARAMN}{Parameter (N)}
#'  \item{AVALC}{Analysis Value (C)}
#'  \item{ANL02FL}{Analysis Flag 02-By Visit Value}
#'  \item{APOBLFL}{Post-Baseline Record Flag}
#'  \item{BASE}{Baseline Value}
#'  \item{BNRIND}{Baseline Reference Range Indicator}
#'  \item{ADTM}{Analysis Date/Time}
#'  \item{CRIT1}{Analysis Criterion 1}
#'  \item{CRIT1FL}{Criterion 1 Evaluation Result Flag}
#'  \item{CRIT2}{Analysis Criterion 2}
#'  \item{CRIT2FL}{Criterion 2 Evaluation Result Flag}
#'  \item{CRIT3}{Analysis Criterion 3}
#'  \item{CRIT3FL}{Criterion 3 Evaluation Result Flag}
#'  \item{ATOXDSCL}{Analysis Toxicity Description Low}
#'  \item{ATOXDSCH}{Analysis Toxicity Description High}
#'  \item{ATOXGRL}{Analysis Toxicity Grade Low}
#'  \item{ATOXGRH}{Analysis Toxicity Grade High}
#'  \item{ATOXGR}{Analysis Toxicity Grade}
#'  \item{ANL06FL}{Analysis Flag 06-Minimum Value}
#'  \item{ANL05FL}{Analysis Flag 05-Worst Tox Grade High}
#'  \item{ANL04FL}{Analysis Flag 04-Worst Value}
#'  \item{ANL03FL}{Analysis Flag 03-Maximum Value}
#'  \item{TRTEMFL}{Treatment Emergent Analysis Flag}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{STUDYID}{Study Identifier}
#'  \item{AGE}{Age}
#'  \item{SEX}{Sex}
#'  \item{RACE_DECODE}{Race}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaefmq}} \code{\link{adaeocmq}} \code{\link{adcm}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets advs
#' @name advs
#' @examples
#' head(data("advs"))
"advs"
