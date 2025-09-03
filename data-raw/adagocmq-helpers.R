#' Derive New Records Based on `ATERMN` Levels
#' 
#' @param df A data frame which can be ADAEOCMQ or ADLB.
#' @param level A single numeric that will be assigned to the derived records.
#' 
#' @return The input data frame with the derived records appended.
#' 
#' @noRd
#' 
#' @name derive_atermn


#' @rdname derive_atermn
#' 
#' For example, if one subject has one record with ATERMN = 121 and another
#' record with ATERMN = 122, and the start date of two those records is within
#' 7 days, then set ATERMN = 12, create a record. If there are multiple
#' combinations satisfy the above criteria, create one record for each
#' combination.
#' 
#' Also applies to ATERMN = 33.
#' 
#' @param levels A numeric vector of two `ATERMN` levels.
#' 
#' @examples
#' df <- dplyr::tibble(
#'   USUBJID = c('a', 'a', 'a'),
#'   ASTDT = as.Date(c("2020-11-01", "2020-11-02", "2020-11-03")),
#'   OCMQNAM = rep("Hypersensitivity", 3),
#'   ATERMN = c(121, 122, 121)
#' )
#' 
#' derive_atermn_1x(df, c(121, 122), 12)
derive_atermn_1x <- function(df, levels, level) {
  ids <- unique(df[["USUBJID"]])
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(ATERMN %in% levels)
    
    if (!all(levels %in% subject[["ATERMN"]])) next
    
    for (i in seq_len(NROW(subject))) {
      current <- subject[i, ]
      
      records_within_7 <- subject[-(1:i), ] |>
        dplyr::filter(abs(ASTDT - current[["ASTDT"]]) <= 7) |>
        dplyr::filter(ATERMN != current[["ATERMN"]]) |>
        dplyr::mutate(ATERMN = level) |>
        dplyr::mutate(ASTDT = min(ASTDT, current[["ASTDT"]]))
      
      df <- dplyr::bind_rows(df, records_within_7)
    }
  }
  
  df
}


#' @rdname derive_atermn
#' 
#' If a subject has more than one records with ATERMN = 231, then set
#' ATERMN = 23, create a record. 
#' 
#' @examples
#' df <- dplyr::tibble(
#'   USUBJID = c("a", "a", "b", "b"),
#'   ATERMN = c(231, 231, 231, NA)
#' )
#' 
#' derive_atermn_23(df)
derive_atermn_23 <- function(df) {
  ids <- unique(df[["USUBJID"]])
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(ATERMN == 231)
    
    if (NROW(subject) <= 1) next
    
    record_23 <- subject[1, ] |>
      dplyr::mutate(ATERMN = 23)
    
    df <- dplyr::bind_rows(df, record_23)
  }
  
  df
}


#' @rdname derive_atermn
#' 
#' If a participant has more than 1 record with ATERMN = 331 and more than 1
#' record with ATERMN = 332, then set ATERMN = 34, create a record. if there
#' are multiple combinations satisfy the above criteria, ONLY keep one
#' combination.
#' 
#' @examples
#' df <- dplyr::tibble(
#'   USUBJID = c(rep("a", 4), rep("b", 3)),
#'   ATERMN = c(331, 331, 332, 332, 331, 331, 332)
#' )
#' 
#' derive_atermn_34(df)
derive_atermn_34 <- function(df) {
  ids <- unique(df[["USUBJID"]])
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(ATERMN %in% c(331, 332)) |>
      dplyr::group_by(USUBJID, ATERMN) |>
      dplyr::mutate(n = dplyr::n())
    
    if (any(subject[["n"]] <= 1)) next
    
    record_34 <- subject |>
      dplyr::ungroup() |>
      dplyr::mutate(n = NULL) |>
      dplyr::slice(1) |>
      dplyr::mutate(ATERMN = 34)
    
    df <- dplyr::bind_rows(df, record_34)
  }
  
  df
}


#' @rdname derive_atermn
#' 
#' if CPK Value >5 x ULN [(ADLB.PARAMCD = 'CK' and ADLB.AVAL> 5*ADLB.ANRHI
#' and ADLB.BASE<=ADLB.ANRHI) and within 3 days, there is no record with
#' (ADLB.PARAM = 'CPK-MB/CPK' and ADLB.AVAL >0.05)], then set ATERMN = 43,
#' create a record
derive_atermn_43 <- function(df) {
  ids <- unique(df[["USUBJID"]])
  records <- df[0, ]
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(PARAMCD == "CK" & AVAL > 5 * ANRHI & BASE <= ANRHI)
    
    if (NROW(subject) == 0) next
    
    dates_within_3 <- subject[["ADT"]]
    dates_within_3 <- c(dates_within_3 + 3, dates_within_3 - 3)
    
    records_within_3 <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(PARAM == "CPK-MB/CPK" & AVAL > 0.05) |>
      dplyr::filter(ADT %in% dates_within_3)
    
    if (NROW(records_within_3) > 0) next
    
    record <- subject[1, ] |>
      dplyr::mutate(ATERMN = 43)
    
    records <- dplyr::bind_rows(records, record)
  }
  
  records
}
