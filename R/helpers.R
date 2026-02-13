#' Check if an argument is NULL or NA
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_null_na <- function(x, name) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) {
    stop(
      sprintf(
        "Given '%s' is NA or NULL, a value must be provided...", name
      )
    )
  }
}

#' Check if argument is a non-missing logical scalar
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_logical <- function(x, name) {
  check_null_na(x, name)
  if (!is.logical(x) || length(x) != 1) {
    stop(
      sprintf(
        "Given '%s' must be a single logical TRUE or FALSE...", name
      )
    )
  }
}

#' Check if argument is a non-missing numeric scalar
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_numeric <- function(x, name) {
  check_null_na(x, name)
  if (!is.numeric(x) || length(x) != 1) {
    stop(
      sprintf(
        "Given '%s' must be a single numeric value...", name
      )
    )
  }
}

#' Check if argument is a non-empty character vector (no NA)
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_character_vector <- function(x, name) {
  check_null_na(x, name)
  if (!is.character(x) || length(x) == 0 || any(is.na(x))) {
    stop(
      sprintf(
        "Given '%s' must be a non-empty character vector without NA...", name
      )
    )
  }
}

#' Check if argument is a numeric vector of exact length
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @param len Expected length of numeric vector
#' @keywords internal
check_numeric_range <- function(x, name, len) {
  check_null_na(x, name)
  if (!is.numeric(x) || length(x) != len || any(is.na(x))) {
    stop(
      sprintf(
        "Given '%s' must be a numeric vector of length %d without NA...",
        name,
        len
      )
    )
  }
}

#' Summary log helper for get_gbif filtering
#'
#' @param log Argument value
#' @param step_name Filtering step name
#' @param before Pre-filtering number of observations
#' @param after Post-filtering number of observations
#' @keywords internal
log_step <- function(log, step_name, before, after) {
  rbind(
    log,
    data.frame(
      step = step_name,
      removed = before - after,
      remaining = after,
      stringsAsFactors = FALSE
    )
  )
}