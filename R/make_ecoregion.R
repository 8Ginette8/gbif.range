### =========================================================================
### make_ecoregion
### =========================================================================
#' Make an ecoregion map based on input environmental variables
#'
#' This function may be used if the in-house ecoregion shapefiles are
#' too coarse for a given geographic region (e.g., for local studies) or
#' a shapefile of finer environmental details is needed. Based on several
#' environmental layers (e.g. climate, soil and land cover), this function
#' can generate a map of environmental regions containing n categories/classes.
#' The classes are calculated with the 'clustering large applications'
#' method (CLARA), which recognize patterns and relationships existing in
#' spatial data, and classify it in clusters.
#'
#' @param env Object of class SpatRaster, RasterBrick or RasterStack of desired
#' resolution and extent defining the study area. Used to generate the map of
#' clusters needed to summarize the environmental space of the study area.
#' @param nclass Numeric, How many number of environmental classes should have
#' the output?
#' @param path Character. Folder path where the output should be saved. Default
#' is none.
#' @param name Character. If 'path' is used, should include the name of the output
#' file (without file extension)
#' @param raster Logical. Whether the output should be a raster layer. Default
#' is FALSE.
#' @param ... Additonnal parameters for the function clara() of the clutser R package.
#' @return A TIFF or SHP file
#' @references
#' Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., &
#' Thuiller, W. (2021). Novel methods to correct for observer and sampling bias
#' in presence‐only species distribution models. Global Ecology and Biogeography,
#' 30(11), 2312-2325.
#' 
#' Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., & Hornik, K. (2021). cluster:
#' Cluster Analysis Basics and Extensions. R package version 2.1.2 — For new features,
#' see the 'Changelog' file (in the package source). https://CRAN.R-project.org/package=cluster
#' 
#' Reynolds, A. P., Richards, G., de la Iglesia, B., & Rayward-Smith, V. J. (2006). Clustering
#' rules: A comparison of partitioning and hierarchical clustering algorithms. Journal of
#' Mathematical Modelling and Algorithms, 5(4), 475–504. doi: 10.1007/s10852-005-9022-1
#' 
#' Schubert, E., & Rousseeuw, P. J. (2019). Faster k-Medoids clustering: Improving the PAM,
#' CLARA, and CLARANS algorithms. In G. Amato, C. Gennaro, V. Oria, & M. Radovanović (Eds.),
#' Similarity search and applications. SISAP 2019. Lecture Notes in Computer Science (Vol.
#' 11807, pp. 171–187). Springer. 
#' 
#' @examples
#' # Load the European Alps extent and a raster of a random resolution
#' data(geo_dat)
#' data(exrst)
#' 
#' # Apply the function by infering 50 classes of environments
#' my.eco = make_ecoregion(rst,50)
#' plot(my.eco)
#' 
#' # Downloading in the European Alps the observations of one plant species
#' obs.arcto = get_gbif("Arctostaphylos alpinus",geo=shp.lonlat)
#' 
#' # Create the range map based on our custom ecoregion
#' range.arcto = get_range("Arctostaphylos alpinus",obs.arcto,my.eco,"EcoRegion",res=20)
#' 
#' # Plot
#' plot(vect(shp.lonlat))
#' plot(range.arcto,add=TRUE,col="darkgreen")
#' points(obs.arcto[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=1)
#' 
#' @export
make_ecoregion=function(env=NULL,nclass=NULL,path="",name="",raster=FALSE,...)
{
    # Check ras input
    if(!(class(env)%in%c("SpatRaster"))) {
      env = rast(env)
    }

    # 'nclass' must be at least equal to '2'
    if (nclass<=1 || is.null(nclass)){
      stop("'nclust' must be equal to 2 or more...!")
    }

    # env must be of a specific class
    if (!(class(env)%in%c("SpatRaster","RasterBrick","RasterStack"))){
      stop("'env' must be an object of class 'SpatRaster','RasterBrick' or 'RasterStack'...!")
    } else if (is.null(env) || nlyr(env)%in%1){
      stop("'ras' must include more than one raster layer...!")
    }

    # First info message
    cat("CLARA algorithm processing...","\n")
  
    # Convert raster in the right CLARA format
    id.toReplace = complete.cases(env[])
    toClara = env[][id.toReplace,]

    # Run CLARA and assign results to right pixels
    blocks = clara(toClara,nclass,...)
    toNew.ras = env[[1]]
    toNew.ras[][id.toReplace] = blocks$clustering
    toNew.ras[][!id.toReplace] = NA

    # Save in raster or shapefile
    if (raster&path=="") {
      return(toNew.ras)

    } else if (raster&path!="") {
      writeRaster(toNew.ras,paste0(path,"/",name,".tif"),overwrite=TRUE,
        datatype="INT4S",gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))
    
    } else {

      # Change raster into a shapefile
      topoly = as.polygons(toNew.ras,dissolve=TRUE)

      # For now, because get_range not yet compatible with Spatvect
      cat("Generating polygons...","\n")
      topoly = as(topoly, "Spatial")
      topoly = suppressWarnings(gBuffer(topoly,byid=TRUE,width=0))
      topoly@data=data.frame(topoly@data,EcoRegion=as.character(1:nrow(topoly@data)))

      if (!raster&path=="") {
        return(topoly)
      
      } else {
        writeVector(topoly,paste0(path,"/",name,".shp"),overwrite=TRUE,
          filetype="ESRI Shapefile",layer=name)
      }
    }
}
