#' @title addisp
#'
#' @description addisp modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 409 rows and 42 variables:
#' \describe{
#'  \item{STUDYID}{Study Identifier}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{SITEID}{Study Site Identifier}
#'  \item{SUBJID}{Subject Identifier for the Study}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{PARAM}{Parameter}
#'  \item{AVALC}{Analysis Value (Character)}
#'  \item{DSSCAT}{Analysis Subcategory}
#'  \item{ASTDT}{Analysis Start Date}
#'  \item{ASTDTC}{Analysis Start Date (Char, ISO8601)}
#'  \item{ASTDY}{Study Day of Analysis Start Date}
#'  \item{AVALCSP}{Specify text for AVALC}
#'  \item{S1EDT}{Date of Last Exposure to Treatment}
#'  \item{S1EDTM}{Datetime of Last Exposure to Treatment}
#'  \item{S1EDY}{S1EDY}
#'  \item{S2EDT}{Date of Last Exposure to Treatment}
#'  \item{S2EDTM}{Datetime of Last Exposure to Treatment}
#'  \item{S2EDY}{S2EDY}
#'  \item{S1SDT}{Date of First Exposure to Treatment}
#'  \item{S1SDTM}{Datetime of First Exposure to Treatment}
#'  \item{S1SDY}{S1SDY}
#'  \item{S2SDT}{Date of First Exposure to Treatment}
#'  \item{S2SDTM}{Datetime of First Exposure to Treatment}
#'  \item{S2SDY}{S2SDY}
#'  \item{PPROTFL}{PPROTFL}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{RANDFL}{RANDFL}
#'  \item{SCRNFL}{SCRNFL}
#'  \item{SCRFFL}{SCRFFL}
#'  \item{RESCRNFL}{RESCRNFL}
#'  \item{FASFL}{FASFL}
#'  \item{ENRLFL}{ENRLFL}
#'  \item{AGE}{Age}
#'  \item{AGEU}{Age Units}
#'  \item{SEX}{Sex}
#'  \item{RACE}{Race}
#'  \item{COUNTRY}{Country}
#'  \item{RFICDT}{RFICDT}
#'  \item{TRT01P}{Planned Treatment for Period 01}
#'  \item{TRT01PN}{TRT01PN}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{TRT01AN}{TRT01AN}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaeocmq}} \code{\link{adagocmq}} \code{\link{adcm}} \code{\link{addili}} \code{\link{addisp}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adpc}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets addisp
#' @name addisp
#' @examples
#'  head(data("addisp"))
"addisp"
