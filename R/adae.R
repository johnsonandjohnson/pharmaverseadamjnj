#' @title adae
#'
#' @description adae modified from pharmaverseadam
#' @source data from pharmaverseadam.
#'
#' @format A data frame with 1191 rows and 83 variables:
#' \describe{
#'  \item{USUBJID}{Unique Subject Identifier}
#'  \item{DOMAIN}{Domain Abbreviation}
#'  \item{AESEQ}{Sequence Number}
#'  \item{AETERM}{Reported Term for the Adverse Event}
#'  \item{AEDECOD}{Dictionary-Derived Term}
#'  \item{AEBODSYS}{Body System or Organ Class}
#'  \item{AEBDSYCD}{Body System or Organ Class Code}
#'  \item{AELLT}{Lowest Level Term}
#'  \item{AELLTCD}{Lowest Level Term Code}
#'  \item{AEPTCD}{Preferred Term Code}
#'  \item{AEHLT}{High Level Term}
#'  \item{AEHLTCD}{High Level Term Code}
#'  \item{AEHLGT}{High Level Group Term}
#'  \item{AEHLGTCD}{High Level Group Term Code}
#'  \item{AESOC}{Primary System Organ Class}
#'  \item{AESOCCD}{Primary System Organ Class Code}
#'  \item{AESTDTC}{Start Date/Time of Adverse Event}
#'  \item{ASTDT}{Analysis Start Date}
#'  \item{ASTDTM}{Analysis Start Date/Time}
#'  \item{ASTDTF}{Analysis Start Date Imputation Flag}
#'  \item{ASTTMF}{Analysis Start Time Imputation Flag}
#'  \item{AEENDTC}{End Date/Time of Adverse Event}
#'  \item{AENDT}{Analysis End Date}
#'  \item{AENDTM}{Analysis End Date/Time}
#'  \item{AENDTF}{Analysis End Date Imputation Flag}
#'  \item{AENTMF}{Analysis End Time Imputation Flag}
#'  \item{ASTDY}{Analysis Start Relative Day}
#'  \item{AESTDY}{Study Day of Start of Adverse Event}
#'  \item{AENDY}{Analysis End Relative Day}
#'  \item{AEENDY}{Study Day of End of Adverse Event}
#'  \item{ADURN}{Analysis Duration (N)}
#'  \item{ADURU}{Analysis Duration Units}
#'  \item{TRTEMFL}{Treatment Emergent Analysis Flag}
#'  \item{AOCCIFL}{1st Max Sev./Int. Occurrence Flag}
#'  \item{AESER}{Serious Event}
#'  \item{AESDTH}{Results in Death}
#'  \item{AESLIFE}{Is Life Threatening}
#'  \item{AESHOSP}{Requires or Prolongs Hospitalization}
#'  \item{AESDISAB}{Persist or Signif Disability/Incapacity}
#'  \item{AESCONG}{Congenital Anomaly or Birth Defect}
#'  \item{AESEV}{Severity/Intensity}
#'  \item{ASEV}{Analysis Severity/Intensity}
#'  \item{ASEVN}{Analysis Severity/Intensity (N)}
#'  \item{AEREL}{Causality}
#'  \item{AREL}{Analysis Causality}
#'  \item{AEACN}{Action Taken with Study Treatment}
#'  \item{AESPID}{Sponsor-Defined Identifier}
#'  \item{AEOUT}{Outcome of Adverse Event}
#'  \item{AESCAN}{Involves Cancer}
#'  \item{AESOD}{Occurred with Overdose}
#'  \item{AEDTC}{Date/Time of Collection}
#'  \item{LDOSEDTM}{End Date/Time of Last Dose}
#'  \item{DOSEON}{Treatment Dose at Record Start}
#'  \item{DOSEU}{Treatment Dose Units}
#'  \item{AETOXGR}{Standard Toxicity Grade}
#'  \item{AETOXGRN}{Standard Toxicity Grade (N)}
#'  \item{DOSEDY}{Day of Study Drug}
#'  \item{AECONTRT}{Concomitant or Additional Trtmnt Given}
#'  \item{CQ01NAM}{Customized Query 01 Name}
#'  \item{CQ02NAM}{Customized Query 02 Name}
#'  \item{CQ03NAM}{Customized Query 03 Name}
#'  \item{AESMIE}{Other Medically Important Serious Event}
#'  \item{ACAT1}{Analysis Category 1}
#'  \item{AOCCFL}{1st Occurrence within Subject Flag}
#'  \item{AOCCPFL}{1st Occurrence within Pref Term Flag}
#'  \item{AOCCSFL}{1st Occurrence of SOC Flag}
#'  \item{AOCTIFL}{1st TE Max Toxicity Grade Flag}
#'  \item{TRT01A}{Actual Treatment for Period 01}
#'  \item{SAFFL}{Safety Population Flag}
#'  \item{AGE}{Age}
#'  \item{SEX}{Sex}
#'  \item{RACE}{Race}
#'  \item{STUDYID}{Study Identifier}
#'  \item{AGEGR1}{Pooled Age Group 1}
#'  \item{TRTEDY}{Treatment Relative End Day}
#'  \item{TRT01P}{Planned Treatment for Period 01}
#'  \item{TRDISCFL}{Treatment Discontinued Flag}
#'  \item{AESHOSPP}{Prolongs Hospitalization}
#'  \item{AESHOSPR}{Requires Hospitalization}
#'  \item{AESCAT}{AESCAT}
#'  \item{SMQ01NAM}{Standardized MedDRA Query 01 Name}
#'  \item{SMQ02NAM}{Standardized MedDRA Query 02 Name}
#'  \item{SMQ03NAM}{Standardized MedDRA Query 03 Name}
#' }
#' @seealso \code{\link{adae}} \code{\link{adae.xpt}} \code{\link{adaeocmq}} \code{\link{adaeocmq.xpt}} \code{\link{adagocmq}} \code{\link{adagocmq.xpt}} \code{\link{adcm}} \code{\link{adcm.xpt}} \code{\link{addili}} \code{\link{addili.xpt}} \code{\link{adeg}} \code{\link{adeg.xpt}} \code{\link{adex}} \code{\link{adex.xpt}} \code{\link{adexsum}} \code{\link{adexsum.xpt}} \code{\link{adlb}} \code{\link{adlb.xpt}} \code{\link{adpc}} \code{\link{adpc.xpt}} \code{\link{adsl}} \code{\link{adsl.xpt}} \code{\link{adttesaf}} \code{\link{adttesaf.xpt}} \code{\link{advs}} \code{\link{advs.xpt}}# nolint
#' @keywords datasets adae
#' @name adae
#' @examples
#'  head(data("adae"))
"adae"
