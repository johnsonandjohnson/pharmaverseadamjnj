#' Derive New Records Based on Rule "ATERMN = 1x"
#' 
#' For example, if one subject has one record with ATERMN = 121 and another
#' record with ATERMN = 122, and the start date of two those records is within
#' 7 days, then set ATERMN = 12, create a record. If there are multiple
#' combinations satisfy the above criteria, create one record for each
#' combination.
#' 
#' Also applies to rule "ATERMN = 33".
#' 
#' @param df A data frame.
#' @param level A single numeric that will be assigned to the derived records.
#' @param levels A numeric vector of two `ATERMN` levels.
#' @returns The input data frame with the derived records appended.
#' 
#' @noRd
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
        dplyr::mutate(ASTDT = min(ASTDT, current[["ASTDT"]])) |>
        dplyr::mutate(
          SRCVALUE = NA,
          SRCVAR = NA,
          SRCDOM = NA,
          SRCSEQ = NA
        )
      
      df <- dplyr::bind_rows(df, records_within_7)
    }
  }
  
  df
}


#' Derive New Records Based on Rule "ATERMN = 23"
#' 
#' If a subject has more than one records with ATERMN = 231, then set
#' ATERMN = 23, create a record. 
#'
#' @param df A data frame.
#' @returns The input data frame with the derived records appended.
#'
#' @noRd
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
    
    record_23 <- subject |>
      dplyr::mutate(ASTDT = min(ASTDT)) |>
      dplyr::slice(1) |>
      dplyr::mutate(ATERMN = 23) |>
      dplyr::mutate(
        SRCVALUE = NA,
        SRCVAR = NA,
        SRCDOM = NA,
        SRCSEQ = NA
      )
    
    df <- dplyr::bind_rows(df, record_23)
  }
  
  df
}


#' Derive New Records Based on Rule "ATERMN = 34"
#' 
#' If a participant has more than 1 record with ATERMN = 331 and more than 1
#' record with ATERMN = 332, then set ATERMN = 34, create a record. if there
#' are multiple combinations satisfy the above criteria, ONLY keep one
#' combination.
#'
#' @param df A data frame.
#' @returns The input data frame with the derived records appended.
#'
#' @noRd
#' @examples
#' df <- dplyr::tibble(
#'   USUBJID = c(rep("a", 4), rep("b", 3), rep("c", 2)),
#'   ATERMN = c(331, 331, 332, 332, 331, 331, 332, 331, 331),
#'   ASTDT = as.Date(c("2023-09-02", rep("2023-09-01", 8)))
#' )
#' 
#' derive_atermn_34(df)
derive_atermn_34 <- function(df) {
  ids <- unique(df[["USUBJID"]])
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(ATERMN %in% c(331, 332))
    
    if (!all(c(331, 332) %in% subject[["ATERMN"]])) next
    
    subject <- subject |>
      dplyr::group_by(ATERMN) |>
      dplyr::mutate(n = dplyr::n())
    
    if (any(subject[["n"]] <= 1)) next
    
    record_34 <- subject |>
      dplyr::ungroup() |>
      dplyr::mutate(ASTDT = min(ASTDT)) |>
      dplyr::mutate(n = NULL) |>
      dplyr::slice(1) |>
      dplyr::mutate(ATERMN = 34) |>
      dplyr::mutate(
        SRCVALUE = NA,
        SRCVAR = NA,
        SRCDOM = NA,
        SRCSEQ = NA
      )
    
    df <- dplyr::bind_rows(df, record_34)
  }
  
  df
}


#' Derive New Records Based on Rule "ATERMN = 44"
#' 
#' If a participant has:
#' 
#' 1. one record from ATERMN = 441 and
#' 2. one record from ATERMN = 442 and
#' 3. one record from ATERMN = 443 and
#' 4. ASTDT for item 1,2, 3 are within 7 days of each other,
#' 
#' then set ATERMN = 44, create a record. If there are multiple combinations
#' satisfy the above criteria, create a record for each combination.
#'
#' @param df A data frame.
#' @returns The input data frame with the derived records appended.
#'
#' @noRd
#' @examples
#' df <- dplyr::tibble(
#'   USUBJID = c(rep("a", 5), "b"),
#'   ASTDT = as.Date(c("2020-11-01", rep("2020-11-02", 5))),
#'   ATERMN = c(441, rep(442, 2), rep(443, 3))
#' )
#' 
#' derive_atermn_44(df)
derive_atermn_44 <- function(df) {
  levels <- c(441, 442, 443)
  ids <- unique(df[["USUBJID"]])
  
  for (id in ids) {
    subject <- df |>
      dplyr::filter(USUBJID == id) |>
      dplyr::filter(ATERMN %in% levels)
    
    if (!all(levels %in% subject[["ATERMN"]])) next
    
    for (i in seq_len(NROW(subject))) {
      current <- subject[i, ]
      current_level <- current[["ATERMN"]]
      current_date <- current[["ASTDT"]]
      
      records_within_7 <- subject[-(1:i), ] |>
        dplyr::filter(abs(ASTDT - current_date) <= 7) |>
        dplyr::filter(ATERMN != current_level)
      
      rest_levels <- setdiff(levels, current_level)
      
      n_combinations <- sum(records_within_7[["ATERMN"]] == rest_levels[1]) *
        sum(records_within_7[["ATERMN"]] == rest_levels[2])
      
      records_within_7 <- records_within_7 |>
        dplyr::mutate(ATERMN = 44) |>
        dplyr::mutate(ASTDT = min(ASTDT, current_date)) |>
        dplyr::slice(rep(1, n_combinations)) |>
        dplyr::mutate(
          SRCVALUE = NA,
          SRCVAR = NA,
          SRCDOM = NA,
          SRCSEQ = NA
        )
      
      df <- dplyr::bind_rows(df, records_within_7)
    }
  }
  
  df
}


#' Create New Records Based on Rule "ATERMN = 43"
#' 
#' If CPK Value > 5 * ULN [(ADLB.PARAMCD = 'CK' and ADLB.AVAL > 5 * ADLB.ANRHI
#' and ADLB.BASE <= ADLB.ANRHI) and within 3 days, there is no record with
#' (ADLB.PARAM = 'CPK-MB/CPK' and ADLB.AVAL > 0.05)], then set ATERMN = 43,
#' create a record.
#' 
#' @param df The ADLB dataset.
#' @returns A data frame that contains only the new records.
#' @noRd
create_atermn_43 <- function(df) {
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
