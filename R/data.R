#' Range area comparison data
#'
#' A dataset containing range area calculations for species comparing 
#' gbif.range-derived polygons and IUCN ranges.
#'
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{species}{Character, the scientific name of the species}
#'   \item{gbif.range.km2}{Numeric, area in square kilometers calculated from gbif.range polygons}
#'   \item{iucn.range.km2}{Numeric, area in square kilometers based on IUCN range data}
#' }
#' @source Calculated using the \code{gbif.range} R package.
#' @docType data
#' @name area_data
NULL