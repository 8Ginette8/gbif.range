### =========================================================================
### obs_filter
### =========================================================================
#' Filter a set of GBIF observations according to a defined grain
#'
#' Whereas the 'grain' parameter in get_gbif() allows GBIF observations to be
#' filtered according to a certain spatial precision, obs_filter() accepts
#' as input a get_gbif() output (one or several species) and filter the
#' observations according to a given grid resolution (one observation
#' per pixel grid kept). This function allows the user to refine the density of
#' GBIF observations according to a defined analysis/study's resolution.
#'
#' @param gbifs One get_gbif() output including one or several species. Note
#' that if GBIF absences are kept in the output(s), the function should be used
#' distinctively for observations and absences.
#' @param grid Object of class SpatRaster, RasterLayer, RasterBrick or
#' RasterStack of desired resolution and extent (WGS84).
#' @return Data frame with two columns named 'x' and 'y' comprising
#' the new set of observations filtered at grid resolution.
#' @example inst/examples/obs_filter_help.R
#' @export
obs_filter=function(gbifs,grid)
{
    # Check 'ras' input
    if(!(class(grid)%in%c("SpatRaster"))) {
      grid = terra::rast(grid)
    }

    # Check number of species in 'get_gbif' output
    n.sp = unique(gbifs$input.search)

    # Loop over species
    out.sp =
    lapply(n.sp,function(x)
    {
      # Extract coordinates of the species
      coords = gbifs[gbifs$input.search%in%x,c("decimalLongitude","decimalLatitude")]

      # Extract related cells
      posP = terra::cellFromXY(grid,as.matrix(coords))

      # Extract one observation per grid cell
      new.oxy = data.frame(Species=x, terra::xyFromCell(grid,unique(posP)))

      # Return
      return(new.oxy)
    })

    # Compile and return
    final.out = do.call("rbind",out.sp)
    return(final.out)
}
