#' @title adttesaf
#'
#' @description adtte_onc modified from pharmaverseadam
#' @source data from adtte_onc from pharmaverseadam
#'
#' @format A data frame with 2032 rows and 9 variables:
#' \describe{
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{PARAM}{Parameter}
#'  \item{AVAL}{Analysis Value}
#'  \item{CNSR}{Censor}
#'  \item{STARTDT}{Start Date}
#'  \item{ADT}{Analysis Date}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{SAFFL}{Safety Population Flag}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaefmq}} \code{\link{adcm}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets adttesaf
#' @name adttesaf
#' @examples
#'  data("adttesaf")
"adttesaf"
