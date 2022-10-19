### =========================================================================
### wsl_obs_filter
### =========================================================================
#' Filter a set of GBIF observations according to a defined grain
#'
#' Filter a set of GBIF observations through a chosen raster grid that defines
#' the resolution of the spatial analysis. Depending on the desired resolution
#' of the analysis, one might want to keep only one GBIF observation per grid
#' cell to avoid e.g., modelling and sampling bias issues.
#'
#' @param wsl.gbif one wsl.gbif output including one or several species
#' @param grid Object of class 'SpatRaster', 'RasterLayer', 'RasterBrick' or
#' 'RasterStack' of desired resolution and extent.
#' @return a data frame with two columns named 'x' and 'y' comprising
#' the new set of observations filtered at grid resolution.
#' @examples
#' 
#' # Load the European Alps extent and a raster of a random resolution
#' data(geo_dat)
#' data(exrst)
#' 
#' # Downloading in the European Alps the observations of Arctostaphylos alpinus
#' obs.arcto = wsl_gbif("Arctostaphylos alpinus",geo=shp.lonlat)
#' plot(shp.lonlat)
#' points(obs.arcto[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=1)
#' 
#' @export
wsl_obs_filter=function(wsl.gbif,grid)
{
    # Check 'ras' input
    if(!(class(grid)%in%c("SpatRaster"))) {
      ras = rast(grid)
    }

    # Check 'wsl.gbif' output
    wsl.gbif$input.search

    # Apply simple filtering
    if (!is.null(a.xy)) {
      
      # Check 'a.xy' input

      if(ncol(a.xy)!=2 || !all(colnames(a.xy)%in%c("x","y"))){
        stop("Supplied points should be a data.frame/matrix with two columns named x and y!")
      }

      # Position presences and absences
      posP=cellFromXY(grid,o.xy)
      posA=cellFromXY(grid,a.xy)

      # Remove absences where we find presences
      posA=posA[!(posA %in% posP)]

      # Extract new presences/absences and regroup
      new.pxy=coordinates(grid)[unique(posP),]
      new.axy=coordinates(grid)[unique(posA),]
      new.oxy=list(new.pxy,new.axy)
      names(new.oxy)=c("Presences","Absences")

    } else {
      # Apply simple filtering
      posCELL=cellFromXY(grid,o.xy)
      new.oxy=coordinates(grid)[unique(posCELL),]
    }

    if (class(new.oxy)[1]%in%"numeric"){
      new.oxy=matrix(new.oxy,ncol=2)
      colnames(new.oxy)=c("x","y")
    }

    return(new.oxy)
}
