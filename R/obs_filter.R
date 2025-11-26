### =========================================================================
### obs_filter
### =========================================================================
#' Filter a set of GBIF observations according to a defined grid resolution
#'
#' The `obs_filter()` function processes a getGBIF object (which can include
#' data for one or multiple species) and filters the observations based on
#' a specified grid resolution. It can retain one observation per grid pixel
#' and/or remove observations from grid pixels that contain fewer than a
#' specified number of records (e.g., fewer than 3 observations per pixel).
#' This function helps users refine the spatial density of GBIF datasets
#' by retaining one observation per grid pixel and/or removing
#' observations from grid pixels with fewer than a user-defined
#' threshold of records.
#'
#' @param gbifs A getGBIF output including one or several species. Note
#' that if GBIF absences are kept in the output(s), the function should be
#' used distinctively for observations and absences.
#' @param grid Object of class SpatRaster, RasterLayer, RasterBrick, or
#' RasterStack of desired resolution and extent (WGS84).
#' @param threshold Optional integer specifying the minimum number of
#' observations per grid pixel to retain. Default is `NULL`, meaning no
#' threshold filtering.
#' @return Data frame with columns 'Species', 'x', and 'y' comprising
#' the new set of observations filtered at grid resolution.
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