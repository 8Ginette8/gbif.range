### ==================================================================
### conv_function (meta)
### ==================================================================
#' Create polygon objects in different bioregions.
#' 
#' Not to be called directly by the user.
#' @param sp_coord spatial coordinates (type "matrix" or "data.frame)
#' @param bwp buffer width parameter (type numeric)
#' @param bipl  number of observation points (type numeric)
#' @param bwpo buffer width parameter for the convex hull (type numeric)
#' @param temp_dir temporary directory name (type character)
#' @param g optional parameter
#' @keywords internal
#' @export
#' @importFrom terra crds buffer aggregate vect crs
#' @importFrom sf st_polygon
conv_function <- function (sp_coord,
                           bwp,
                           bipl,
                           bwpo,
                           temp_dir,
                           g=NULL){
  
  # Preps and convert degrees to meters
  x = terra::crds(sp_coord)
  row.names(x) = 1:nrow(x)
  bwp_bipl_m = (bwp+(nrow(x)-1)*bipl)*111111
  bwpo_m = bwpo*111111

  # Check number of observations points in each bioregion, if <3 create point buffer
  if (nrow(data.frame(x))<3) { 
    rtn = terra::buffer(terra::aggregate(sp_coord),width=bwp_bipl_m)
    return(rtn) 
    
  } else {
    
    # Test if points are on line, if yes create point buffer around points
    is_line = 0
    for (i in 2:(nrow(x)-1))
    {
      dxc = x[i,1]-x[i-1,1]
      dyc = x[i,2]-x[i-1,2]
      dx1 = x[i+1,1]-x[i-1,1]
      dy1 = x[i+1,2]-x[i-1,2]  
      is_line[i-1] = dxc*dy1-dyc*dx1
    }
    
    if (all(abs(is_line) == 0)) { 
      
      # Print
      cat('bioreg=',g,nrow(x),'points laying on one line. Using buffer width of ',bwp_bipl_m/1000,'km','\n')

      # Buffer
      rtn = terra::buffer(terra::aggregate(sp_coord),width=bwp_bipl_m)
      
      # Out
      return(rtn)

    } else { 

      if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive=TRUE)
      }

      # If they are not on a line, create a convex hull      
      arg_file = paste0('QJ Fx TO ', file.path(temp_dir, 'vert.txt'))
      vert0 = geometry::convhulln(x, arg_file)
      vert1 = scan(file.path(temp_dir,'vert.txt'),quiet=T)
      file.remove(file.path(temp_dir,'vert.txt'))
      vert2 = (vert1+1)[-1]
      FE_vert = row.names(x)[vert2]
      coord_conv = x[FE_vert,]

      # Polygonizing
      coord_conv = rbind(coord_conv,coord_conv[1,])
      P1 = sf::st_polygon(list(coord_conv))
      P1 = terra::vect(P1)
      terra::crs(P1) = terra::crs(sp_coord)
      rtn = terra::buffer(P1,width=bwpo_m)
      terra::crs(rtn) = terra::crs(sp_coord) # Safety measure...

      # Out
      return(rtn)
    }
  }
} 
