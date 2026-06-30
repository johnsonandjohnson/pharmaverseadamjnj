#' @title adishum
#'
#' @description adishum modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 3856 rows and 56 variables:
#' \describe{
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{ISSEQ}{Sequence Number}
#'  \item{ISBDAGNT}{Binding Agent}
#'  \item{ISSTRESC}{Character Result/Finding in Std Format}
#'  \item{ISSTRESN}{Numeric Results/Findings in Std. Units}
#'  \item{ISNAM}{Vendor Name}
#'  \item{VISITNUM}{Visit Number}
#'  \item{VISIT}{Visit Name}
#'  \item{ISDTC}{Date/Time of Collection}
#'  \item{STUDYID}{Study Identifier}
#'  \item{SUBJID}{Subject Identifier for the Study}
#'  \item{SITEID}{Study Site Identifier}
#'  \item{AGE}{Age}
#'  \item{SEX}{Sex}
#'  \item{RACE}{Race}
#'  \item{IMFL}{Immunogenicity Population Flag}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{TRTSDT}{Date of First Exposure to Treatment}
#'  \item{TRTSDTM}{Datetime of First Exposure to Treatment}
#'  \item{TRT01P}{Planned Treatment for Period 01}
#'  \item{TRT01PN}{Planned Treatment for Period 01 (N)}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{TRT01AN}{Actual Treatment for Period 01 (N)}
#'  \item{ADTM}{Analysis Datetime}
#'  \item{TRTA}{TRTA}
#'  \item{PARQUAL}{Parameter Qualifier}
#'  \item{ISTESTCD}{Immunogenicity Test/Exam Short Name}
#'  \item{ISTEST}{Immunogenicity Test or Examination Name}
#'  \item{ISSPEC}{Specimen Type}
#'  \item{ISCAT}{Category for Immunogenicity Test}
#'  \item{ISSCAT}{ISSCAT}
#'  \item{ISTSTDTL}{ISTSTDTL}
#'  \item{PARCAT1}{PARCAT1}
#'  \item{PARAMCON}{PARAMCON}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{PARAM}{Parameter Description}
#'  \item{AVAL}{Analysis Value}
#'  \item{AVALC}{Analysis Value (C)}
#'  \item{AVALCAT1}{Analysis Value Category 1}
#'  \item{IMEVFL}{IS Evaluable Population Flag}
#'  \item{AVISITN}{Analysis Visit (N)}
#'  \item{AVISIT}{Analysis Visit}
#'  \item{VISITDY}{Planned Study Day of Visit}
#'  \item{ADT}{Analysis Date}
#'  \item{ADY}{Analysis Relative Day}
#'  \item{ABLFL}{Baseline Record Flag}
#'  \item{ANL01FL}{Analysis Flag 01}
#'  \item{ANL02FL}{Analysis Flag 02}
#'  \item{ANL03FL}{Analysis Flag 03-Inf SR}
#'  \item{ANL07FL}{Analysis Flag 07-Inj SR}
#'  \item{ANL04FL}{Analysis Flag 04-Severe Inf SR}
#'  \item{ANL05FL}{Analysis Flag 05-Serious Inf SR}
#'  \item{ANL06FL}{Analysis Flag 06-Infus SR - DC}
#'  \item{ANL08FL}{Analysis Flag 08-Severe Inj SR}
#'  \item{ANL09FL}{Analysis Flag 09-Serious Inj SR}
#'  \item{ANL10FL}{Analysis Flag 10-Inject SR - DC}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaecomp}} \code{\link{adaeocmq}} \code{\link{adagocmq}} \code{\link{adcm}} \code{\link{addili}} \code{\link{addisp}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adishum}} \code{\link{adlb}} \code{\link{adpc}} \code{\link{adsl}} \code{\link{adslcomp}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets adishum
#' @name adishum
#' @examples
#'  head(data("adishum"))
"adishum"
