### =========================================================================
### wsl.obs.filter
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
#' # Load my binary observations species data
#' data(var_select_XYtest)
#' data(exrst)
#' 
#' # wsl.obs.filter(): example for the first species
#'    # Loading observations: presences and absences
#' presences = coordinates(mySP[[1]])[myPA[[1]] %in% "1",]
#' absences = coordinates(mySP[[1]])[myPA[[1]] %in% "0",]
#'
#'    # Loading grid
#' 
#' r.layer = rst[[1]]
#'
#'    # To filter observations by the grid only
#' 
#' pres.filtered = wsl.obs.filter(presences,grid=r.layer)
#' abs.filtered = wsl.obs.filter(absences,grid=r.layer)
#'
#'    # To filter observations by the grid & remove abs. in cells where we also find pres.
#' 
#' PresAbs.filtered = wsl.obs.filter(presences,absences,r.layer)
#'
#'    # Count presences (same filtering)
#' 
#' nrow(PresAbs.filtered[[1]])
#' nrow(pres.filtered)
#'
#'    # Count Absences (filtering plus removal of duplicated absences)
#' 
#' nrow(abs.filtered)
#' nrow(PresAbs.filtered[[2]])
#'
#'    # Visual
#' 
#' par(mfrow=c(1,2))
#' plot(presences)
#' plot(pres.filtered)
#'
#' par(mfrow=c(1,2))
#' plot(absences)
#' plot(abs.filtered)
#' 
#' @export
wsl.obs.filter=function(wsl.gbif,grid)
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
