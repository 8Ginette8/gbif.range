#' Range Area Comparison Data
#'
#' A dataset containing range-area estimates for species, comparing
#' \code{gbif.range}-derived polygons with corresponding IUCN ranges.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{species}{Character string with the scientific name of the species.}
#'   \item{gbif.range.km2}{Numeric area in square kilometers calculated from
#'   \code{gbif.range} polygons.}
#'   \item{iucn.range.km2}{Numeric area in square kilometers calculated from
#'   IUCN range data.}
#' }
#' @source Calculated with the \code{gbif.range} package.
#' @docType data
#' @name area_data
NULL
