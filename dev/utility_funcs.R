#' Find package datasets that contain specified variable names
#'
#' @description
#' Search all documented datasets in a package (the list returned by
#' `data(package = pkg)`) and report which datasets contain one or more
#' variable names you provide.
#'
#' @param pkg Character(1). Name of the installed package to search (e.g. "pharmaverseadam").
#' @param vars Character vector. One or more variable names to look for (e.g. c("PARAMCD")).
#' @param match_all Logical(1). If TRUE, only return datasets that contain all variables
#'   in `vars`. If FALSE (default), return datasets that contain any of the variables.
#'
#' @return A tibble with columns:
#'   - dataset: dataset name in the package
#'   - matched_vars: list-column of matched variable names found in that dataset
#'   - data: list-column with the full dataset object (or NULL if load_data = FALSE)
#'
#' @examples
#' \dontrun{
#' # Find datasets that contain both PARAMCD and PARAM
#' find_datasets_with_vars("pharmaverseadam", c("PARAMCD", "PARAM"), match_all = TRUE)
#' }
#'
find_datasets_with_vars <- function(pkg, vars, match_all = FALSE, load_data = TRUE) {
  stopifnot(is.character(pkg), length(pkg) == 1, is.character(vars), length(vars) >= 1)
  ds <- utils::data(package = pkg)$results[, "Item"]

  out_dataset <- character()
  out_matched <- list()
  out_data <- list()

  for (d in ds) {
    obj <- tryCatch(
      utils::getFromNamespace(d, pkg),
      error = function(e1) {
        tryCatch(
          utils::getExportedValue(pkg, d),
          error = function(e2) {
            tmp_env <- new.env()
            tryCatch(
              {
                utils::data(list = d, package = pkg, envir = tmp_env)
                get(d, envir = tmp_env)
              },
              error = function(e3) NULL
            )
          }
        )
      }
    )

    if (is.null(obj)) {
      next
    }

    obj_names <- names(obj) # works for both data.frame and list; NULL for atomic
    if (is.null(obj_names)) {
      next
    }

    matched <- intersect(vars, obj_names)
    keep <- if (match_all) length(matched) == length(vars) else length(matched) > 0
    if (!keep) {
      next
    }

    out_dataset <- c(out_dataset, d)
    out_matched[[length(out_matched) + 1]] <- matched
    out_data[[length(out_data) + 1]] <- if (load_data) obj else NULL
  }

  if (length(out_dataset) == 0) {
    return(list(
      all_datasets = list(),
      summary = tibble::tibble(dataset_name = character(), matched_vars = character())
    ))
  }

  names(out_data) <- out_dataset

  list(
    all_datasets = out_data,
    summary = tibble::tibble(
      dataset_name = out_dataset,
      matched_vars = vapply(out_matched, paste, character(1), collapse = ", ")
    )
  )
}

#' Convert RDA to SAS7BDAT
#'
#' Loads all data frames from an .rda file and writes each as a .sas7bdat file.
#'
#' @param rda_path Path to the .rda file.
#' @param output_dir Output directory. Defaults to the same directory as the input file.
#'
#' @return Invisibly returns the output file paths.
#' @examples
#' convert_rda_to_sas7bdat("path/to/your/file.rda")
#' convert_rda_to_sas7bdat("path/to/your/file.rda", output_dir = "path/to/output")
convert_rda_to_sas7bdat <- function(rda_path, output_dir = dirname(rda_path)) {
  env <- new.env()
  load(rda_path, envir = env)

  objs <- ls(env)
  out_paths <- character(0)
  for (obj_name in objs) {
    df <- get(obj_name, envir = env)
    if (is.data.frame(df)) {
      df <- df |>
        dplyr::mutate(
          dplyr::across(
            where(is.factor),
            as.character
          )
        )
      out_path <- file.path(
        output_dir,
        paste0(obj_name, ".sas7bdat")
      )
      haven::write_sas(df, out_path)
      message("Written: ", out_path)
      out_paths <- c(out_paths, out_path)
    }
  }
  invisible(out_paths)
}

#' Convert XPT to SAS7BDAT
#'
#' Reads an .xpt file and writes it as a .sas7bdat file.
#'
#' @param xpt_path Path to the .xpt file.
#' @param output_dir Output directory. Defaults to the same directory as the input file.
#'
#' @return Invisibly returns the output file path.
#' @examples
#' convert_xpt_to_sas7bdat("path/to/your/file.xpt")
#' convert_xpt_to_sas7bdat("path/to/your/file.xpt", output_dir = "path/to/output")
convert_xpt_to_sas7bdat <- function(xpt_path, output_dir = dirname(xpt_path)) {
  df <- haven::read_xpt(xpt_path)
  obj_name <- tools::file_path_sans_ext(basename(xpt_path))
  out_path <- file.path(output_dir, paste0(obj_name, ".sas7bdat"))
  haven::write_sas(df, out_path)
  message("Written: ", out_path)
  invisible(out_path)
}
