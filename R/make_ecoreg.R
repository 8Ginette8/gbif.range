### =========================================================================
### make_ecoreg
### =========================================================================
#' Build Custom Ecoregions from Environmental Layers
#'
#' Cluster multi-layer environmental data to create a custom ecoregion map that
#' can be used directly in \code{get_range()}.
#'
#' @param env Environmental raster stack. Accepted classes are
#' \code{SpatRaster}, \code{RasterBrick}, and \code{RasterStack}.
#' @param nclass Numeric number of environmental classes to create.
#' @param path Optional directory where the output should be written. Leave
#' empty to return the result directly.
#' @param name Output file name without extension when \code{path} is used.
#' @param format Output format for polygon results: \code{"SpatVector"}
#' (default) or \code{"sf"}.
#' @param raster Logical. Should the output remain a raster instead of being
#' converted to polygons? Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{cluster::clara()}.
#' @details This function is useful when the packaged ecoregion layers are too
#' coarse for a study area or when a custom environmental regionalization is
#' needed. Clusters are created with the CLARA algorithm on the multivariate
#' environmental space represented by \code{env}.
#' @return If \code{path == ""}, returns the generated raster or polygon object.
#' Otherwise, writes the output to disk as a GeoTIFF or Shapefile.
#' @references
#' Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., &
#' Thuiller, W. (2021). Novel methods to correct for observer and sampling bias
#' in presence‐only species distribution models. Global Ecology and
#' Biogeography, 30(11), 2312-2325.
#' 
#' Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., & Hornik, K. (2021).
#' cluster: Cluster Analysis Basics and Extensions. R package version 2.1.2 —
#' For new features, see the 'Changelog' file (in the package source).
#' https://CRAN.R-project.org/package=cluster
#' 
#' Reynolds, A. P., Richards, G., de la Iglesia, B., & Rayward-Smith, V. J.
#' (2006). Clustering rules: A comparison of partitioning and hierarchical
#' clustering algorithms. Journal of Mathematical Modelling and Algorithms,
#' 5(4), 475–504. doi: 10.1007/s10852-005-9022-1
#' 
#' Schubert, E., & Rousseeuw, P. J. (2019). Faster k-Medoids clustering:
#' Improving the PAM, CLARA, and CLARANS algorithms. In G. Amato, C. Gennaro,
#' V. Oria, & M. Radovanović (Eds.), Similarity search and applications.
#' SISAP 2019. Lecture Notes in Computer Science (Vol. 11807, pp. 171–187).
#' Springer. 
#' @example inst/examples/make_ecoreg_help.R
#' @importFrom terra rast nlyr writeRaster as.polygons vect buffer writeVector
#' @importFrom sf st_make_valid st_as_sf
#' @importFrom cluster clara
#' @importFrom stats complete.cases
#' @export
make_ecoreg <- function(env = NULL,
                          nclass = NULL,
                          path = "",
                          name = "",
                          format = "SpatVector",
                          raster = FALSE,
                          ...)
{
    ######################################################
    ### Stop messages
    ######################################################


    # General
    check_numeric(nclass, "nclass")
    check_character_vector(path, "path")
    check_character_vector(name, "name")
    check_logical(raster, "raster")

    # env must be of a specific class
    if (!(class(env) %in% c("SpatRaster", "RasterBrick", "RasterStack"))){
      stop(
        paste(
          "'env' must be an object of class 'SpatRaster',",
          "'RasterBrick' or 'RasterStack'...!"
        )
      )
      
    } else if (is.null(env) || terra::nlyr(env) %in%1 ){
      stop("'ras' must include more than one raster layer...!")
    }

    # Check ras input
    if(!(class(env) %in% c("SpatRaster"))) {
      env <- terra::rast(env)
    }

    # 'nclass' must be at least equal to '2'
    if (nclass <= 1 || is.null(nclass)){
      stop("'nclust' must be equal to 2 or more...!")
    }


    ######################################################
    ### Make ecoregions
    ######################################################


    # First info message
    cat("CLARA algorithm processing...","\n")
  
    # Convert raster in the right CLARA format
    id.toReplace <- stats::complete.cases(env[])
    toClara <- env[][id.toReplace, ]

    # Run CLARA and assign results to right pixels
    blocks <- cluster::clara(toClara, nclass, ...)
    toNew.ras <- env[[1]]
    toNew.ras[][id.toReplace] <- blocks$clustering
    toNew.ras[][!id.toReplace] <- NA

    # Save in raster or shapefile
    if (raster&path == "") {
      return(toNew.ras)

    } else if (raster&path != "") {
      terra::writeRaster(
        x = toNew.ras,
        filename = paste0(path ,"/", name, ".tif"),
        overwrite = TRUE,
        datatype = "INT4S",
        gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
      )

    } else {

      # Change raster into a shapefile
      topoly <- terra::as.polygons(toNew.ras, dissolve = TRUE)

      # For now, because get_range not yet compatible with Spatvect
      cat("Generating polygons...","\n")
      names(topoly) <- "CLARA"
      topoly$EcoRegion <- as.character(seq_len(nrow(data.frame(topoly))))

      # Testing if polygons are valid and correct if not
      topoly <- terra::makeValid(topoly)
      topoly <- terra::buffer(topoly, width = 0) 
      
      if (format == "sf") {
        # Return sf object
        topoly <- sf::st_as_sf(topoly)
        if (!raster & path == "") {
          return(topoly)
        } else {
          sf::st_write(
            obj = topoly,
            dsn = file.path(path, paste0(name, ".shp")),
            layer = name,
            driver = "ESRI Shapefile",
            delete_layer = TRUE
            )
          }
      } else {
        if (!raster & path == "") {
          return(topoly)
        } else {
          terra::writeVector(
            x = topoly,
            filename = file.path(path, paste0(name, ".shp")),
            overwrite = TRUE,
            filetype = "ESRI Shapefile",
            layer = name
          )
        }
      }
    }
}
