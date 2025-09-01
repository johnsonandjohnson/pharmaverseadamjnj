#' Derive New Records Based on `ATERMN` Levels
#' 
#' @param df A data frame which can be ADAEOCMQ or ADLB.
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
#' @param levels A numeric vector of two `ATERMN` levels.
#' @param level A single numeric that will be assigned to the derived records.
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
