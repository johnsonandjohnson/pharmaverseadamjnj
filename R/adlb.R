#' @title adlb
#'
#' @description adlb modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 83640 rows and 156 variables:
#' \describe{
#'  \item{STUDYID}{Study Identifier}
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{SUBJID}{Subject Identifier for the Study}
#'  \item{SITEID}{Study Site Identifier}
#'  \item{COUNTRY}{Country}
#'  \item{DOMAIN}{Domain Abbreviation}
#'  \item{RFSTDTC}{Subject Reference Start Date/Time}
#'  \item{RFENDTC}{Subject Reference End Date/Time}
#'  \item{RFXSTDTC}{Date/Time of First Study Treatment}
#'  \item{RFXENDTC}{Date/Time of Last Study Treatment}
#'  \item{RFPENDTC}{Date/Time of End of Participation}
#'  \item{SCRFDT}{Screen Failure Date}
#'  \item{FRVDT}{Final Retrieval Visit Date}
#'  \item{DTHDTC}{Date/Time of Death}
#'  \item{DTHADY}{Relative Day of Death}
#'  \item{DTHFL}{Subject Death Flag}
#'  \item{LDDTHELD}{Elapsed Days from Last Dose to Death}
#'  \item{LDDTHGR1}{Last Dose to Death - Days Elapsed Grp 1}
#'  \item{DTH30FL}{Death Within 30 Days of Last Trt Flag}
#'  \item{DTHA30FL}{Death After 30 Days from Last Trt Flag}
#'  \item{DTHDOM}{Domain for Date of Death Collection}
#'  \item{DTHB30FL}{Death Within 30 Days of First Trt Flag}
#'  \item{ASEQ}{Analysis Sequence Number}
#'  \item{REGION1}{Geographic Region 1}
#'  \item{DMDTC}{Date/Time of Collection}
#'  \item{DMDY}{Study Day of Collection}
#'  \item{AGE}{Age}
#'  \item{AGEU}{Age Units}
#'  \item{SEX}{Sex}
#'  \item{RACE}{Race}
#'  \item{RACEGR1}{Pooled Race Group 1}
#'  \item{ETHNIC}{Ethnicity}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{ARM}{Description of Planned Arm}
#'  \item{ARMCD}{Planned Arm Code}
#'  \item{ACTARM}{Description of Actual Arm}
#'  \item{ACTARMCD}{Actual Arm Code}
#'  \item{TRTP}{Planned Treatment}
#'  \item{TRTA}{Actual Treatment}
#'  \item{TRT01P}{Planned Treatment for Period 01}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{TRTSDT}{Date of First Exposure to Treatment}
#'  \item{TRTSDTM}{Datetime of First Exposure to Treatment}
#'  \item{TRTSTMF}{Time of First Exposure Imput. Flag}
#'  \item{TRTEDT}{Date of Last Exposure to Treatment}
#'  \item{TRTEDTM}{Datetime of Last Exposure to Treatment}
#'  \item{TRTETMF}{Time of Last Exposure Imput. Flag}
#'  \item{EOSSTT}{End of Study Status}
#'  \item{EOSDT}{End of Study Date}
#'  \item{RFICDTC}{Date/Time of Informed Consent}
#'  \item{RANDDT}{Date of Randomization}
#'  \item{LSTALVDT}{Date Last Known Alive}
#'  \item{TRTDURD}{Total Treatment Duration (Days)}
#'  \item{DTHDT}{Date of Death}
#'  \item{DTHDTF}{Date of Death Imputation Flag}
#'  \item{DTHCAUS}{Cause of Death}
#'  \item{DTHCGR1}{Cause of Death Reason 1}
#'  \item{ADT}{Analysis Date}
#'  \item{ADY}{Analysis Relative Day}
#'  \item{AVISIT}{Analysis Visit}
#'  \item{AVISITN}{Analysis Visit (N)}
#'  \item{PARAM}{Parameter}
#'  \item{PARAMCD}{Parameter Code}
#'  \item{PARAMN}{Parameter (N)}
#'  \item{PARCAT1}{Parameter Category 1}
#'  \item{AVAL}{Analysis Value}
#'  \item{AVALC}{Analysis Value (C)}
#'  \item{BASE}{Baseline Value}
#'  \item{BASEC}{Baseline Value (C)}
#'  \item{BASETYPE}{Baseline Type}
#'  \item{CHG}{Change from Baseline}
#'  \item{PCHG}{Percent Change from Baseline}
#'  \item{R2BASE}{Ratio to Baseline}
#'  \item{R2ANRLO}{Ratio of Analysis Val compared to ANRLO}
#'  \item{R2ANRHI}{Ratio of Analysis Val compared to ANRHI}
#'  \item{SHIFT1}{Shift from Baseline to Analysis Value}
#'  \item{SHIFT2}{Shift from Baseline to Overall Grade}
#'  \item{DTYPE}{Derivation Type}
#'  \item{ATOXGR}{Analysis Toxicity Grade}
#'  \item{BTOXGR}{Baseline Toxicity Grade}
#'  \item{ANRIND}{Analysis Reference Range Indicator}
#'  \item{BNRIND}{Baseline Reference Range Indicator}
#'  \item{ANRLO}{Analysis Normal Range Lower Limit}
#'  \item{ANRHI}{Analysis Normal Range Upper Limit}
#'  \item{ATOXGRL}{Analysis Toxicity Grade Low}
#'  \item{ATOXGRH}{Analysis Toxicity Grade High}
#'  \item{BTOXGRL}{Baseline Toxicity Grade Low}
#'  \item{BTOXGRH}{Baseline Toxicity Grade High}
#'  \item{ATOXDSCL}{Analysis Toxicity Description Low}
#'  \item{ATOXDSCH}{Analysis Toxicity Description High}
#'  \item{ABLFL}{Baseline Record Flag}
#'  \item{ANL01FL}{Analysis Flag 01}
#'  \item{ONTRTFL}{On Treatment Record Flag}
#'  \item{LVOTFL}{Last Value On Treatment Record Flag}
#'  \item{LBSEQ}{Sequence Number}
#'  \item{LBTESTCD}{Lab Test or Examination Short Name}
#'  \item{LBTEST}{Lab Test or Examination Name}
#'  \item{LBCAT}{Category for Lab Test}
#'  \item{LBORRES}{Result or Finding in Original Units}
#'  \item{LBORRESU}{Original Units}
#'  \item{LBORNRLO}{Reference Range Lower Limit in Orig Unit}
#'  \item{LBORNRHI}{Reference Range Upper Limit in Orig Unit}
#'  \item{LBSTRESC}{Character Result/Finding in Std Format}
#'  \item{LBSTRESN}{Numeric Result/Finding in Standard Units}
#'  \item{LBSTRESU}{Standard Units}
#'  \item{LBSTNRLO}{Reference Range Lower Limit-Std Units}
#'  \item{LBSTNRHI}{Reference Range Upper Limit-Std Units}
#'  \item{LBNRIND}{Reference Range Indicator}
#'  \item{LBBLFL}{Baseline Flag}
#'  \item{VISITNUM}{Visit Number}
#'  \item{VISIT}{Visit Name}
#'  \item{VISITDY}{Planned Study Day of Visit}
#'  \item{LBDTC}{Date/Time of Specimen Collection}
#'  \item{LBDY}{Study Day of Specimen Collection}
#'  \item{TRT01PN}{Planned Treatment for Period 01 (N)}
#'  \item{TRT01AN}{Actual Treatment for Period 01 (N)}
#'  \item{AVALU}{Analysis Value - Units}
#'  \item{ANL02FL}{Analysis Record Flag 02-Analysis Value}
#'  \item{TRTEMFL}{Treatment Emergent Analysis Flag}
#'  \item{COUNTRY_DECODE}{Country}
#'  \item{RACE_DECODE}{Race Description}
#'  \item{ETHNIC_DECODE}{Ethnicity Description}
#'  \item{PARCAT2}{Parameter Category 2}
#'  \item{PARCAT3}{Parameter Category 3}
#'  \item{PARCAT4}{Parameter Category 4}
#'  \item{PARCAT5}{Parameter Category 5}
#'  \item{PARCAT6}{Parameter Category 6}
#'  \item{MCRIT2ML}{Multi-Response Criterion 2 Evaluation}
#'  \item{MCRIT1ML}{Multi-Response Criterion 1 Evaluation}
#'  \item{MCRIT1MN}{Multi-Response Criterion 1 Eval (N)}
#'  \item{MCRIT2MN}{Multi-Response Criterion 2 Eval (N)}
#'  \item{MCRIT1}{Analysis Multi-Response Criterion 1}
#'  \item{MCRIT2}{Analysis Multi-Response Criterion 2}
#'  \item{APOBLFL}{Post-Baseline Record Flag}
#'  \item{LBSTNRHQ}{Reference Limit Higher}
#'  \item{LBSTNRLQ}{Reference Limit Lower}
#'  \item{ATOXGRN}{Analysis Toxicity Grade (Numeric)}
#'  \item{ADTM}{Analysis Date/Time}
#'  \item{ATPT}{Analysis Timepoint}
#'  \item{LBCLSIG}{Clinically Significant}
#'  \item{TR01SDT}{Start Date of Treatment for Period 01}
#'  \item{TR01EDT}{End Date of Treatment for Period 01}
#'  \item{LBSPEC}{Specimen Type}
#'  \item{LBFAST}{Fasting Status}
#'  \item{LBNAM}{Laboratory Name}
#'  \item{ANL03FL}{Analysis Record Flag 03 - Protocol Visit}
#'  \item{ANL04FL}{Analysis Flag 04}
#'  \item{ANL05FL}{Analysis Flag 05}
#'  \item{ANL06FL}{Analysis Flag 06}
#'  \item{ANL07FL}{Analysis Flag 07}
#'  \item{ANL08FL}{Analysis Flag 08}
#'  \item{ANL09FL}{Analysis Flag 09}
#'  \item{ANL10FL}{Analysis Flag 10}
#'  \item{ANL14FL}{Analysis Flag 14}
#'  \item{ANL15FL}{Analysis Flag 15}
#'  \item{ANL16FL}{Analysis Flag 16}
#' }
#' @seealso \code{\link{adae}} \code{\link{adaeocmq}} \code{\link{adagocmq}} \code{\link{adcm}} \code{\link{addili}} \code{\link{adeg}} \code{\link{adex}} \code{\link{adexsum}} \code{\link{adlb}} \code{\link{adpc}} \code{\link{adsl}} \code{\link{adttesaf}} \code{\link{advs}}# nolint
#' @keywords datasets adlb
#' @name adlb
#' @examples
#' head(data("adlb"))
"adlb"
