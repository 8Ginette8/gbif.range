### =========================================================================
### Meta function: create polygon objects in different bioregions
### =========================================================================
#' Create polygon objects in different bioregions
#'
#' Not to be called directly by the user
#' @export
conv_function <- function (sp_coord,
                           bwp,
                           bipl,
                           bwpo,
                           temp_dir,
                           g=NULL){
  
  x=sp_coord@coords
  if (is.null(row.names(x))) {row.names(x)=1:nrow(x)}

  if (nrow(x)<3){ #check number of observations points in each bioregion, if <3 create point buffer
    rtn=gBuffer(sp_coord,width=bwp+(nrow(x)-1)*bipl)@polygons[[1]]
    return(rtn) 
    
  } else {
    
    #test if points are on line, if yes create point buffer around points
    is_line <- 0
    for(i in 2:(nrow(x)-1)){
      dxc <- x[i,1]-x[i-1,1]
      dyc <- x[i,2]-x[i-1,2]
      dx1 <- x[i+1,1]-x[i-1,1]
      dy1 <- x[i+1,2]-x[i-1,2]  
      is_line[i-1] <- dxc*dy1-dyc*dx1
    }
    
    if(all(abs(is_line)==0)){ 
      cat('Bioreg=', g, nrow(x), 'points laying on one line. Using buffer width of ', bwp+(nrow(x)-1)*bipl, '\n')

      rtn=gBuffer(sp_coord,width=bwp+(nrow(x)-1)*bipl)@polygons[[1]]
      return(rtn)
    }
    else { #if they are not on a line, create a convex hull
      
      if(!dir.exists(temp_dir)){
        dir.create(temp_dir, recursive=T)
      }
      
      arg_file <- paste0('QJ Fx TO ', file.path(temp_dir, 'vert.txt'))
      vert0<-convhulln(x, arg_file)
      vert1<-scan(file.path(temp_dir,'vert.txt'),quiet=T)
      file.remove(file.path(temp_dir,'vert.txt'))
      vert2<-(vert1+1)[-1]
      FE_vert<-row.names(x)[vert2]
      coord_conv <- x[FE_vert,]  
      
      coord_conv <-rbind(coord_conv,coord_conv[1,])
      P1 <- Polygons(srl=list(Polygon(coord_conv,hole=FALSE)),ID="PolygA")
      P1 <- SpatialPolygons(Srl=list(P1),proj4string=sp_coord@proj4string)
      
      return(gBuffer(P1,width=bwpo)@polygons[[1]])
    }
  }
} 
