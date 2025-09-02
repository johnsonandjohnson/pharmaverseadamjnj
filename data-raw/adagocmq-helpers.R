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
#' @param A single character that specifies the criterion.
#' @return The derived records.
#' 
#' @examples
#' df <- dplyr::tibble(
#'   OCMQNAM = rep("Hypersensitivity", 2),
#'   HYPSCAT = c("A", "B")
#' )
#' 
#' derive_atermn(df, 'OCMQNAM == "Hypersensitivity" & HYPSCAT == "A"', 11)
derive_atermn <- function(df, rule, level) {
  df |>
    dplyr::filter(eval(parse(text = rule))) |>
    dplyr::mutate(ATERMN = level)
}


#' @rdname derive_atermn
#' 
#' For example, if one subject has one record with ATERMN = 121 and another
#' record with ATERMN = 122, and the start date of two those records is within
#' 7 days, then set ATERMN = 12, create a record. If there are multiple
#' combinations satisfy the above criteria, create one record for each
#' combination.
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
derive_atermn_24 <- function(df) {
  # Some variables should include or exclude some keywords
  cmindc_include <- paste(c(
    "diab",
    "mellitus",
    "hyperglyc",
    "glucose",
    "dibet",
    "dieb"
  ), collapse = "|")
  
  cmindc_exclude <- paste(c(
    "prophyla",
    "prevent",
    "insipidus",
    "hyperglycerid",
    "low blood glucose",
    "low glucose",
    "low blood sugar",
    "low sugar",
    "low afternoon blood glucose",
    "low morning blood glucose"
  ), collapse = "|")
  
  cmclas_include <- paste(c(
    "gliptin",
    "glutide",
    "diabet",
    "glitaz",
    "glucose lowering",
    "glucosidas",
    "dipeptidyl",
    "sulfonyl",
    "DPP",
    "guanide",
    "GLP",
    "glucagon-like",
    "metform",
    "gliflozin",
    "insulin",
    "sodium-glucose",
    "SGLT",
    "thiazolid"
  ), collapse = "|")
  
  cmclas_exclude <- "sex hormone"
  
  # Derive `ATERMN`
  df |>
    dplyr::filter(
      ASTDT >= TRTSDT &
        grepl(cmindc_include, CMINDC, ignore.case = TRUE) &
        !grepl(cmindc_exclude, CMINDC, ignore.case = TRUE) &
        grepl(cmclas_include, CMCLAS, ignore.case = TRUE) &
        !grepl(cmclas_exclude, CMCLAS, ignore.case = TRUE)
    ) |>
    dplyr::mutate(ATERMN = 24)
}
