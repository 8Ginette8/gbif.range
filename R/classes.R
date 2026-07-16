### =========================================================================
### Set get_range class
### =========================================================================
#' Reference Class for \code{get_range()} Results
#'
#' Stores the original arguments used to build a range map and the resulting
#' spatial output.
#'
#' @return A generator object of reference class \code{"getRange"}, used to
#' instantiate objects with two fields: \code{init.args} (the original
#' arguments used to build the range map) and \code{rangeOutput} (the
#' resulting spatial output).
#' @export
getRange <- setRefClass("getRange",
                        fields = list(
                          init.args = "ANY",
                          rangeOutput = "ANY"
                        )
)

### =========================================================================
### Set get_gbif class
### =========================================================================
#' Constructor for \code{getGBIF} Objects
#'
#' Wrap a data frame in the \code{getGBIF} class used by downstream functions
#' in the package.
#'
#' @param df A data frame containing GBIF occurrence records.
#' @return An object of class \code{c("getGBIF", "data.frame")}.
#' @export
getGBIF <- function(df) {
    if (!is.data.frame(df)) {
        stop("Input must be a data.frame")
    }
    structure(df, class = c("getGBIF", "data.frame"))
}
