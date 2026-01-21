### =========================================================================
### Set get_range class
### =========================================================================
#' Class to store get_range outuput
#'
#' @keywords internal
getRange <- setRefClass("getRange",
                        fields = list(
                          init.args = "ANY",
                          rangeOutput = "ANY"
                        )
)

### =========================================================================
### Set get_gbif class
### =========================================================================
#' Class to store get_gbif outuput
#'
#' @keywords internal
getGBIF <- function(df) {
    if (!is.data.frame(df)) {
        stop("Input must be a data.frame")
    }
    structure(df, class = c("getGBIF", "data.frame"))
}