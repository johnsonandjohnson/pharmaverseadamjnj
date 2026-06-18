#' @title addili
#'
#' @description addili modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 10476 rows and 71 variables:
#' \describe{
#'  \item{STUDYID}{Study Identifier}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{PARAM}{Parameter}
#'  \item{ADT}{Analysis Date}
#'  \item{ADTM}{Analysis Date/Time}
#'  \item{ADY}{Analysis Relative Day}
#'  \item{AVISIT}{Analysis Visit}
#'  \item{AVISITN}{Analysis Visit (N)}
#'  \item{ATPT}{Analysis Timepoint}
#'  \item{TRTP}{Planned Treatment}
#'  \item{TRTA}{Actual Treatment}
#'  \item{AVAL}{Analysis Value}
#'  \item{BASE}{Baseline Value}
#'  \item{ANRIND}{Analysis Reference Range Indicator}
#'  \item{BNRIND}{Baseline Reference Range Indicator}
#'  \item{ANRLO}{Analysis Normal Range Lower Limit}
#'  \item{ANRHI}{Analysis Normal Range Upper Limit}
#'  \item{ABLFL}{Baseline Record Flag}
#'  \item{APOBLFL}{Post-Baseline Record Flag}
#'  \item{ONTRTFL}{On Treatment Record Flag}
#'  \item{LBSEQ}{Sequence Number}
#'  \item{VISIT}{Visit Name}
#'  \item{VISITNUM}{Visit Number}
#'  \item{LBNAM}{Laboratory Name}
#'  \item{DTYPE}{Derivation Type}
#'  \item{TRT01P}{Planned Treatment for Period 01}
#'  \item{TRT01PN}{Planned Treatment for Period 01 (N)}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{TRT01AN}{Actual Treatment for Period 01 (N)}
#'  \item{TRTSDT}{Date of First Exposure to Treatment}
#'  \item{TRTSDTM}{Datetime of First Exposure to Treatment}
#'  \item{TRTEDT}{Date of Last Exposure to Treatment}
#'  \item{TRTEDTM}{Datetime of Last Exposure to Treatment}
#'  \item{AGE}{Age}
#'  \item{AGEU}{Age Units}
#'  \item{AGEGR1}{Pooled Age Group 1}
#'  \item{AGEGR1N}{Pooled Age Group 1 (N)}
#'  \item{SEX}{Sex}
#'  \item{RACE}{Race}
#'  \item{COUNTRY}{Country}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{FASFL}{Full Analysis Set Population Flag}
#'  \item{RANDFL}{Randomized Population Flag}
#'  \item{DTHFL}{Subject Death Flag}
#'  \item{SITEID}{Study Site Identifier}
#'  \item{SUBJID}{Subject Identifier for the Study}
#'  \item{PARCAT1}{Parameter Category 1}
#'  \item{R2ANRHI}{Ratio to ANR Upper Limit}
#'  \item{CRIT1}{Analysis Criterion 1}
#'  \item{CRIT1FL}{Criterion 1 Evaluation Result Flag}
#'  \item{CRIT2}{Analysis Criterion 2}
#'  \item{CRIT2FL}{Criterion 2 Evaluation Result Flag}
#'  \item{CRIT3}{Analysis Criterion 3}
#'  \item{CRIT3FL}{Criterion 3 Evaluation Result Flag}
#'  \item{CRIT4}{Analysis Criterion 4}
#'  \item{CRIT4FL}{Criterion 4 Evaluation Result Flag}
#'  \item{CRIT5}{Analysis Criterion 5}
#'  \item{CRIT5FL}{Criterion 5 Evaluation Result Flag}
#'  \item{ANL04FL}{Analysis 04 Record Flag}
#'  \item{ANL06FL}{Analysis 06 Record Flag}
#'  \item{ANL05FL}{Analysis 05 Record Flag}
#'  \item{ANL07FL}{Analysis 07 Record Flag}
#'  \item{ANL02FL}{Analysis 02 Record Flag}
#'  \item{ANL03FL}{Analysis 03 Record Flag}
#'  \item{AVALC}{Analysis Value (C)}
#'  \item{AVALCAT1}{Analysis Value Category 1}
#'  \item{AVALCA1N}{Analysis Value Category 1 (N)}
#'  \item{AVALCAT2}{Analysis Value Category 2}
#'  \item{AVALCA2N}{Analysis Value Category 2 (N)}
#'  \item{ASEQ}{Analysis Sequence Number}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaecomp}} \code{\link{adaeocmq}} \code{\link{adagocmq}} \code{\link{adcm}} \code{\link{addili}} \code{\link{addisp}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adishum}} \code{\link{adlb}} \code{\link{adpc}} \code{\link{adsl}} \code{\link{adslcomp}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets addili
#' @name addili
#' @examples
#'  head(data("addili"))
"addili"
