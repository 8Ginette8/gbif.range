### ==================================================================
### conv_function (meta)
### ==================================================================
#' Low-Level Polygon Builder Used by \code{get_range()}
#' 
#' Internal helper that creates buffered polygons from clustered occurrence
#' points inside a single ecoregion.
#'
#' @param sp_coord Spatial coordinates for one cluster of occurrences.
#' @param bwp Numeric buffer width, in degrees, used for isolated points.
#' @param bipl Numeric increment, in degrees, added to the point buffer for
#' linear clusters.
#' @param bwpo Numeric buffer width, in degrees, applied to convex hull
#' polygons.
#' @param temp_dir Character string giving the temporary directory used for
#' intermediate convex-hull files.
#' @param g Optional ecoregion identifier used in status messages.
#' @keywords internal
#' @importFrom terra crds buffer aggregate vect crs
#' @importFrom sf st_polygon
#' @export
conv_function <- function (sp_coord,
                           bwp,
                           bipl,
                           bwpo,
                           temp_dir,
                           g = NULL){
  
  # Preps and convert degrees to meters
  x <- terra::crds(sp_coord)
  row.names(x) <- seq_len(nrow(x))
  bwp.bipl.m <- (bwp + (nrow(x)-1) * bipl) * 111111
  bwpo.m <- bwpo * 111111

  # Number of observations points in each ecoregion, if <3 create point buffer
  if (nrow(data.frame(x)) < 3) { 
    rtn <- terra::buffer(
      x = terra::aggregate(sp_coord),
      width = bwp.bipl.m
    )
    return(rtn) 
    
  } else {
    
    # Test if points are on line, if yes create point buffer around points
    is.line <- 0
    for (i in 2:(nrow(x)-1))
    {
      dxc <- x[i, 1] - x[i-1, 1]
      dyc <- x[i, 2] - x[i-1, 2]
      dx1 <- x[i+1, 1] - x[i-1, 1]
      dy1 <- x[i+1, 2] - x[i-1, 2]  
      is.line[i-1] <- dxc * dy1 - dyc * dx1
    }
    
    if (all(abs(is.line) == 0)) { 
      
      # Print
      cat('\n','ecoreg=',g,nrow(x),
        'points lying on one line. Using buffer width of ',
        bwp.bipl.m/1000,'km','\n')

      # Buffer
      rtn <- terra::buffer(
        x = terra::aggregate(sp_coord),
        width = bwp.bipl.m
      )
      
      # Out
      return(rtn)

    } else { 

      if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
      }

      # If they are not on a line, create a convex hull      
      arg.file <- paste0('QJ Fx TO ', file.path(temp_dir, 'vert.txt'))
      vert0 <- geometry::convhulln(x, arg.file)
      vert1 <- scan(file.path(temp_dir, 'vert.txt'),quiet = TRUE)
      file.remove(file.path(temp_dir, 'vert.txt'))
      vert2 <- (vert1 + 1)[-1]
      fe.vert <- row.names(x)[vert2]
      coord.conv <- x[fe.vert, ]

      # Polygonizing
      coord.conv <- rbind(coord.conv, coord.conv[1, ])
      P1 <- sf::st_polygon(list(coord.conv))
      P1 <- terra::vect(P1)
      terra::crs(P1) <- terra::crs(sp_coord)
      rtn <- terra::buffer(
        x = P1,
        width = bwpo.m
      )
      terra::crs(rtn) <- terra::crs(sp_coord) # Safety measure...

      # Out
      return(rtn)
    }
  }
} 
