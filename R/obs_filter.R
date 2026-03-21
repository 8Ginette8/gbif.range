### =========================================================================
### obs_filter
### =========================================================================
#' Filter GBIF Records by Grid Cell
#'
#' Reduce the spatial density of a \code{getGBIF} object by retaining a single
#' record per grid cell and, optionally, removing cells with too few records.
#'
#' @param gbifs A \code{getGBIF} object containing one or more species.
#' @param grid Raster defining the target spatial resolution and extent.
#' Accepted classes are \code{SpatRaster}, \code{RasterLayer},
#' \code{RasterBrick}, and \code{RasterStack}.
#' @param threshold Optional integer specifying the minimum number of records a
#' cell must contain to be retained.
#' @details The function first collapses each species to one occurrence per grid
#' cell. If \code{threshold} is supplied, cells with fewer than that many
#' original records are then discarded.
#' @return A data frame with the columns \code{Species}, \code{x}, and
#' \code{y}, representing the filtered coordinates.
#' @example inst/examples/obs_filter_help.R
#' @importFrom terra rast cellFromXY xyFromCell
#' @export
obs_filter <- function(gbifs, grid, threshold = NULL) {


  #####################################################
  ### Stop messages
  #####################################################


  # gbifs
  if (!methods::is(gbifs,"getGBIF")) {
    stop("Given 'gbifs' must be a getGBIF object.")
  }

  # grid
  if (!class(grid) %in% c("SpatRaster", "RasterLayer", "RasterBrick")){
    stop("Given 'grid' must be a 'SpatRaster', 'RasterLayer' or 'RasterBrick'")
  }


  #####################################################
  ### Filter observations
  #####################################################


  # Check 'grid' input
  if (!inherits(grid, "SpatRaster")) {
    grid <- terra::rast(grid)
  }

  # Check number of species in the getGBIF object
  n.sp <- unique(gbifs$input_search)

  # Loop over species
  out.sp <- lapply(n.sp, function(x) {
    
    # Extract coordinates of the species
    coords <- gbifs[gbifs$input_search %in% x,
          c("decimalLongitude", "decimalLatitude")]

    # Extract related cells
    posP <- terra::cellFromXY(grid, as.matrix(coords))

    # Extract one observation per grid cell
    unique_posP <- unique(posP)
    new_oxy <- data.frame(
      Species = x,
      terra::xyFromCell(grid, unique_posP)
    )

    # Apply threshold filtering if specified
    if (!is.null(threshold)) {
      cell_counts <- table(posP)
      filtered_posP <- unique_posP[unique_posP %in%
                names(cell_counts[cell_counts >= threshold])]
      terxy <- terra::cellFromXY(grid, new_oxy[, c("x", "y")])
      new_oxy <- new_oxy[terxy %in% filtered_posP, ]
    }

    # Return
    return(new_oxy)
  })

  # Compile and return
  final.out <- do.call("rbind", out.sp)
  return(final.out)
}
