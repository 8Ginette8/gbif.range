#' Merge the polygons of a species range map
#'
#' @description
#' \code{get_range()} returns one polygon per occupied ecoregion. Those
#' polygons do not overlap, so their areas are additive and the ecoregion
#' origin of every piece stays visible. \code{merge_range()} dissolves them
#' into a single polygon, or into groups given by \code{by}. The dissolve is
#' kept out of \code{get_range()} because collapsing the
#' ecoregions discards the \code{ecoreg_name} attribute. Rasterized output
#' (\code{format = "SpatRaster"}) never needs this function.
#'
#' @param x A \code{getRange} object returned by \code{get_range()}, or the
#' \code{SpatVector} or \code{sf} object held in its \code{rangeOutput} field.
#' @param by Character. Optional name of a field in \code{x} used to group the
#' polygons: one feature is returned per distinct value. Note that on the
#' default output of \code{get_range()} this is already the case, so
#' \code{by = ecoreg_name} returns the input unchanged. It is useful for
#' coarser groupings, such as a biome field, or for objects whose polygons were
#' not built by \code{get_range()}. When \code{NULL} (default), every polygon is
#' dissolved into a single feature.
#' @param format Character. Output format. One of \code{"SpatVector"}
#' (default) or \code{"sf"}.
#' @return A \code{SpatVector} or \code{sf} object holding the dissolved range.
#' With \code{by} supplied, the grouping field is retained along with any other
#' field that is constant within each group. With \code{by = NULL} the result
#' carries no attributes, since no attribute of the original polygons describes
#' the single merged geometry.
#' @seealso \code{\link{get_range}}() to build the range map.
#' @example inst/examples/merge_range_help.R
#' @importFrom methods is
#' @importFrom terra vect aggregate geom crs
#' @export
merge_range <- function(x,
                        by = NULL,
                        format = c("SpatVector", "sf")) {

  format <- match.arg(format)

  # Accept the getRange container or the geometry directly
  if (methods::is(x, "getRange")) {
    x <- x$rangeOutput
  }

  if (methods::is(x, "SpatRaster")) {
    stop(
      paste(
        "'x' is already rasterized and holds no polygons to merge.",
        "Re-run get_range() with format = 'SpatVector' or 'sf'."
      )
    )
  }

  if (methods::is(x, "sf")) {
    x <- terra::vect(x)
  }

  if (!methods::is(x, "SpatVector")) {
    stop("'x' must be a getRange object, a SpatVector, or an sf object...")
  }

  if (nrow(x) == 0) {
    stop("'x' contains no polygons...")
  }

  # Grouping field
  if (!is.null(by)) {
    if (!methods::is(by, "character") || length(by) != 1) {
      stop("'by' must be a single character string...")
    }
    if (!by %in% names(x)) {
      stop(paste0("'", by, "' is not a field of 'x'..."))
    }
    out <- terra::aggregate(x, by = by, count = FALSE)
  } else {
    # Attributes are meaningless once every polygon is dissolved into one, so
    # rebuild from geometry alone. This also sidesteps terra issue #1710, where
    # a SpatVector returned by erase() rejects the column aggregate() would add.
    out <- terra::aggregate(terra::vect(terra::geom(x), type = "polygons"))
    terra::crs(out) <- terra::crs(x)
  }

  if (format == "sf") {
    out <- sf::st_as_sf(out)
  }

  return(out)
}
