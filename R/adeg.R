#' @title adeg
#'
#' @description adeg modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 11844 rows and 70 variables:
#' \describe{
#'  \item{DOMAIN}{Domain Abbreviation}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{EGSEQ}{Sequence Number}
#'  \item{EGTESTCD}{ECG Test or Examination Short Name}
#'  \item{EGTEST}{ECG Test or Examination Name}
#'  \item{EGORRES}{Result or Finding in Original Units}
#'  \item{EGORRESU}{Original Units}
#'  \item{EGSTRESC}{Character Result/Finding in Std Format}
#'  \item{EGSTRESN}{Numeric Result/Finding in Standard Units}
#'  \item{EGSTRESU}{Standard Units}
#'  \item{EGSTAT}{Completion Status}
#'  \item{EGLOC}{Lead Location Used for Measurement}
#'  \item{EGBLFL}{Baseline Flag}
#'  \item{VISITNUM}{Visit Number}
#'  \item{VISIT}{Visit Name}
#'  \item{VISITDY}{Planned Study Day of Visit}
#'  \item{EGDTC}{Date/Time of ECG}
#'  \item{EGDY}{Study Day of ECG}
#'  \item{EGTPT}{Planned Time Point Name}
#'  \item{EGTPTNUM}{Planned Time Point Number}
#'  \item{EGELTM}{Planned Elapsed Time from Time Point Ref}
#'  \item{EGTPTREF}{Time Point Reference}
#'  \item{ADTM}{Analysis Datetime}
#'  \item{ATMF}{Analysis Time Imputation Flag}
#'  \item{ADY}{Analysis Relative Day}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{AVAL}{Analysis Value}
#'  \item{AVALC}{Analysis Value (C)}
#'  \item{ADT}{Analysis Date}
#'  \item{ATPTN}{Analysis Timepoint (N)}
#'  \item{ATPT}{Analysis Timepoint}
#'  \item{AVISIT}{Analysis Visit}
#'  \item{AVISITN}{Analysis Visit (N)}
#'  \item{DTYPE}{Derivation Type}
#'  \item{ONTRTFL}{On Treatment Record Flag}
#'  \item{ANRLO}{Analysis Normal Range Lower Limit}
#'  \item{ANRHI}{Analysis Normal Range Upper Limit}
#'  \item{ANRIND}{Analysis Reference Range Indicator}
#'  \item{BASETYPE}{Baseline Type}
#'  \item{ABLFL}{Baseline Record Flag}
#'  \item{BASEC}{Baseline Value (C)}
#'  \item{CHG}{Change from Baseline}
#'  \item{PCHG}{Percent Change from Baseline}
#'  \item{ANL01FL}{Analysis Flag 01-Analysis Value}
#'  \item{TRTP}{Planned Treatment}
#'  \item{TRTA}{Actual Treatment}
#'  \item{ASEQ}{Analysis Sequence Number}
#'  \item{AVALCAT1}{Analysis Value Category 1}
#'  \item{AVALCA1N}{Analysis Value Category 1 (N)}
#'  \item{CHGCAT1}{Change from Baseline Category 1}
#'  \item{CHGCAT1N}{Change from Baseline Category 1 (N)}
#'  \item{PARAM}{Parameter}
#'  \item{PARAMN}{Parameter (N)}
#'  \item{TRTEMFL}{Treatment Emergent Analysis Flag}
#'  \item{ANL02FL}{Analysis Flag 02-By Visit Value}
#'  \item{ANL03FL}{Analysis Flag 03-Maximum Value}
#'  \item{APOBLFL}{Post-Baseline Record Flag}
#'  \item{CRIT1}{Analysis Criterion 1}
#'  \item{CRIT1FL}{Criterion 1 Evaluation Result Flag}
#'  \item{CRIT2}{Analysis Criterion 2}
#'  \item{CRIT2FL}{Criterion 2 Evaluation Result Flag}
#'  \item{BASE}{Baseline Value}
#'  \item{BNRIND}{Baseline Reference Range Indicator}
#'  \item{BASECAT1}{Baseline Category 1}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{STUDYID}{Study Identifier}
#'  \item{AGE}{Age}
#'  \item{SEX}{Sex}
#'  \item{RACE_DECODE}{Race}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaefmq}} \code{\link{adcm}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets adeg
#' @name adeg
#' @examples
#'  head(data("adeg"))
"adeg"
