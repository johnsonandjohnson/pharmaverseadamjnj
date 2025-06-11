#' @title adsl
#'
#' @description adsl modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 306 rows and 106 variables:
#' \describe{
#'  \item{STUDYID}{Study Identifier}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{SUBJID}{Subject Identifier for the Study}
#'  \item{RFSTDTC}{Subject Reference Start Date/Time}
#'  \item{RFENDTC}{Subject Reference End Date/Time}
#'  \item{RFXSTDTC}{Date/Time of First Study Treatment}
#'  \item{RFXENDTC}{Date/Time of Last Study Treatment}
#'  \item{RFICDTC}{Date/Time of Informed Consent}
#'  \item{RFPENDTC}{Date/Time of End of Participation}
#'  \item{DTHDTC}{Date/Time of Death}
#'  \item{DTHFL}{Subject Death Flag}
#'  \item{SITEID}{Study Site Identifier}
#'  \item{AGE}{Age}
#'  \item{AGEU}{Age Units}
#'  \item{SEX}{Sex}
#'  \item{RACE}{Race}
#'  \item{ETHNIC}{Ethnicity}
#'  \item{ARMCD}{Planned Arm Code}
#'  \item{ARM}{Description of Planned Arm}
#'  \item{ACTARMCD}{Actual Arm Code}
#'  \item{ACTARM}{Description of Actual Arm}
#'  \item{COUNTRY}{Country}
#'  \item{DMDTC}{Date/Time of Collection}
#'  \item{DMDY}{Study Day of Collection}
#'  \item{TRT01P}{Planned Treatment for Period 01}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{TRTSDTM}{Datetime of First Exposure to Treatment}
#'  \item{TRTSTMF}{Time of First Exposure Imput. Flag}
#'  \item{TRTEDTM}{Datetime of Last Exposure to Treatment}
#'  \item{TRTETMF}{Time of Last Exposure Imput. Flag}
#'  \item{TRTSDT}{Date of First Exposure to Treatment}
#'  \item{TRTEDT}{Date of Last Exposure to Treatment}
#'  \item{TRTDURD}{Total Treatment Duration (Days)}
#'  \item{SCRFDT}{Screen Failure Date}
#'  \item{EOSDT}{End of Study Date}
#'  \item{EOSSTT}{End of Study Status}
#'  \item{FRVDT}{Final Retrievel Visit Date}
#'  \item{RANDDT}{Date of Randomization}
#'  \item{DTHDT}{Date of Death}
#'  \item{DTHDTF}{Date of Death Imputation Flag}
#'  \item{DTHADY}{Relative Day of Death}
#'  \item{LDDTHELD}{Elapsed Days from Last Dose to Death}
#'  \item{DTHCAUS}{DTHCAUS}
#'  \item{DTHDOM}{DTHDOM}
#'  \item{DTHCGR1}{DTHCGR1}
#'  \item{LSTALVDT}{Date Last Known Alive}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{RACEGR1}{Pooled Race Group 1}
#'  \item{AGEGR1}{Pooled Age Group 1}
#'  \item{REGION1}{Geographic Region 1}
#'  \item{LDDTHGR1}{Last Dose to Death - Days Elapsed Grp 1}
#'  \item{DTH30FL}{Death Within 30 Days of Last Trt Flag}
#'  \item{DTHA30FL}{Death After 30 Days from Last Trt Flag}
#'  \item{DTHB30FL}{Death Within 30 Days of First Trt Flag}
#'  \item{TRT01PN}{Planned Treatment for Period 01 (N)}
#'  \item{TRT01AN}{Actual Treatment for Period 01 (N)}
#'  \item{AGEGR1N}{Pooled Age Group 1 (N)}
#'  \item{SEX_DECODE}{Sex}
#'  \item{WEIGHTBL}{Weight (kg)}
#'  \item{WGTGR1N}{Weight Group 1 (N)}
#'  \item{WGTGR1}{Weight Group 1}
#'  \item{HEIGHTBL}{Height (cm)}
#'  \item{BSABL}{Body surface area (m2)}
#'  \item{BMIBL}{Body mass index (kg/m2)}
#'  \item{BMIBLG1N}{BMI at Baseline Group 1 (N)}
#'  \item{BMIBLG1}{BMI at Baseline Group 1}
#'  \item{COUNTRY_DECODE}{Country}
#'  \item{RACE_DECODE}{Race}
#'  \item{RFICDT}{Date of Informed Consent}
#'  \item{ETHNIC_DECODE}{Ethnicity}
#'  \item{STRAT1R}{Strat Factor 1 Value Used for Rand}
#'  \item{STRAT2R}{Strat Factor 2 Value Used for Rand}
#'  \item{RANUM}{Randomization Number}
#'  \item{RANDDTM}{Datetime of Randomization}
#'  \item{EOTSTT}{End of Treatment Status}
#'  \item{DCTREAS}{Reason for Discontinuation of Treatment}
#'  \item{LTVISIT}{Last Treatment Visit}
#'  \item{DCTREASP}{Reason Specify for Discont of Treatment}
#'  \item{DCTDT}{End of Study Date}
#'  \item{DCSREAS}{Reason for Discontinuation from Study}
#'  \item{DCSREASP}{Reason Spec for Discont from Study}
#'  \item{LSVISIT}{Last Study Visit}
#'  \item{TRTEDY}{Treatment Relative End Day}
#'  \item{SCRNFL}{Screened Population Flag}
#'  \item{SCRFFL}{Screen Failure Flag}
#'  \item{DCSCREEN}{Reason for Discont During Screening}
#'  \item{ENRLFL}{Enrolled Population Flag}
#'  \item{RANDFL}{Randomized Flag}
#'  \item{ITTFL}{Intent-To-Treat Population Flag}
#'  \item{FASFL}{Full Analysis Set Population Flag}
#'  \item{PPROTFL}{Per-Protocol Population Flag}
#'  \item{LSTSVDT}{Last Subject Visit (SV) Date}
#'  \item{EOSDY}{Study Day of Study Termination}
#'  \item{UNBLNDFL}{Subject Blind Broken}
#'  \item{RESCRNFL}{Re-screened Flag}
#'  \item{DTHTRTFL}{Death on Treatment Flag}
#'  \item{DTHCAUSP}{Cause Spec for Death}
#'  \item{DTHAFTFL}{Death After 30 Days of Last Treatment}
#'  \item{DTH60TFL}{Death Within 60 Days of First Treatment}
#'  \item{UNBLNDDY}{Study Day of Unblinding}
#'  \item{UNBREAS}{Reason For Unblinding}
#'  \item{LDOSE}{Last Dose}
#'  \item{LDOSU}{Last Dose Unit}
#'  \item{DTHTERM}{Reported Cause of Death}
#'  \item{LDSTODTH}{Days from Last Dose to Death}
#'  \item{DTHDY}{Study Day of Death}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaefmq}} \code{\link{adcm}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets adsl
#' @name adsl
#' @examples
#' \dontrun{
#' data("adsl")
#' }
"adsl"
