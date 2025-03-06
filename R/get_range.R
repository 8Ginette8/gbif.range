### ==================================================================
### get_range
### ==================================================================
#' Create a species range map based on a get_gbif() output
#' 
#' Estimates species ranges based on occurrence data (GBIF or not) and ecoregions
#' (related functions 'make_ecoregion'). It first deletes outliers from the observation
#' dataset and then creates a polygon (convex hull) with a user specified buffer around
#' all the observations of one ecoregion. If there is only one observation in an ecoregion,
#' a buffer around this point will be created. If all points in an ecoregion are on a line,
#' the function will also create a buffer around these points, however, the buffer size
#' increases with the number of points in the line. Finally, also note that in case of
#' too many records, get_range can be used with a sub-sample of species observations to
#' ensure a faster polygon process and/or to overcome potential RAM crash of the function.
#' 
#' @param occ_coord a get_gbif() output or a data.frame containing two columns named
#' "decimalLongitude" and "decimalLatitude".
#' @param bioreg  'SpatialPolygonsDataFrame', 'SpatVector' or 'sf' object containg different
#' ecoregions (convex hulls will be classified on a bioreg basis) and of CRS WGS84. Note
#' that this parameter may be fed with an external, generated (function make_ecoregion) or
#' in-house ecoregion shapefile. Four shapefiles can be downloaded with the library (see functions read_bioreg()
#' and others): 'eco_terra' (for terrestrial species; Nature conservancy version adapted from Olson & al. 2001),
#' 'eco_marine' and 'eco_hd_marine' (for marine species; Spalding & al. 2007, 2012) and 'eco_fresh' (for freshwater
#' species; Abell & al. 2008). For marine species, 'eco_terra' may also be used if the user wants to represent
#' the terrestrial range of species that also partially settle on mainland. For fresh water species, same may be
#' done if the user considers that terrestrial ecoregions should be more representative of the species ecology.
#' @param bioreg_name Character. One ecoregion level/category name from the 'bioreg' polygon must be supplied.
#' E.g., very detailed level of 'eco_terra' is 'ECO_NAME'. Note that default applies if a make_ecoregion()
#' polygon is provided as 'bioreg'.
#' @param degrees_outlier Numeric. Distance threshold (degrees) for outlier classification.
#' If the nearest minimal distance to the next point is larger than this threshold, it will be
#' considered as an outlier.
#' @param clustered_points_outlier Numeric. Maximum number of points which are closer to each other
#' than the degrees_outlier, but should still be considered as outliers.
#' @param buffer_width_point Numeric. Buffer (in degrees) which will be applied around single
#' observations.
#' @param buffer_increment_point_line Numeric. How much should the buffer be increased for each point
#' on a line.
#' @param buffer_width_polygon Numeric. Buffer (in degrees) which will be applied around distribution
#' polygons (for each ecoregion).
#' @param dir_temp Character. Where should the temporary text file for the convex hull be saved?
#' (text file will be deleted again). Default value is \code{tempdir()}.
#' @param raster Logical. Should the output be a unified raster? Default is TRUE
#' @param res Numeric. If raster = TRUE, which resolution of the output in degrees (1° = ~111 km at
#' the equator). Default is 0.01 (~1.1 km). It is important to note that the highest achievable resolution
#' of the output will depend on its 'bioreg' precision, e.g., a species range output can reach the same
#' resolution of the rasters used to create a 'make_ecoregion' object.
#' @param verbose Logical. Should progession be printed?
#' @details Ecoregions cover relatively large areas of land or water, and contain characteristic,
#' geographically distinct assemblages of natural communities sharing a large majority of species,
#' dynamics, and environmental conditions. The biodiversity of flora, fauna and ecosystems that
#' characterise an ecoregion tends to be distinct from that of other ecoregions
#' (https://en.wikipedia.org/wiki/Ecoregion).
#' 
#' Each ecoregion shapefile has one or more categories, which describe more or less precisely the
#' ecoregion world distribution (from the more to the less detailed):
#' 
#' - eco_terra has three different levels: 'ECO_NAME', 'WWF_MHTNAM' and 'WWF_REALM2'.
#' - eco_fresh has only one: 'ECOREGION'.
#' - eco_marine and eco_hd_marine (very coastal-precise version) contains three distinct levels:
#' 'ECOREGION', 'PROVINCE' and 'REALM'.
#' 
#' @return A 'SpatVector' or 'SpatRaster' + the list of all arguments used in the function.
#' @references
#' Oskar Hagen, Lisa Vaterlaus, Camille Albouy, Andrew Brown, Flurin Leugger, Renske E. Onstein,
#' Charles Novaes de Santana, Christopher R. Scotese, Loïc Pellissier. (2019) Mountain building,
#' climate cooling and the richness of cold-adapted plants in the Northern Hemisphere. Journal of
#' Biogeography. doi: 10.1111/jbi.13653
#' 
#' Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J. S., ... &
#' Pellissier, L. (2022). An integrated high resolution mapping shows congruent biodiversity
#' patterns of Fagales and Pinales. New Phytologist, 235(2), 759-772 10.1111/nph.18158
#' 
#' Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N.,
#' Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J.,
#' Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem,
#' K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth.
#' Bioscience 51(11):933-938. doi: 10.1641/0006-3568(2001)051
#' 
#' The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types, Biogeographical Realms
#' and The Nature Conservancy Terrestrial Assessment Units. GIS layers developed by The Nature
#' Conservancy with multiple partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986.
#' Cambridge (UK): The Nature Conservancy.
#' 
#' Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A. Ferdaña, Max
#' Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana, Sara A. Lourie, Kirsten D.
#' Martin, Edmund McManus, Jennifer Molnar, Cheri A. Recchia, James Robertson, Marine
#' Ecoregions of the World: A Bioregionalization of Coastal and Shelf Areas, BioScience,
#' Volume 57, Issue 7, July 2007, Pages 573–583. doi: 10.1641/B570707
#' 
#' Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012). Pelagic provinces of the world:
#' a biogeographic classification of the world’s surface pelagic waters. Ocean & Coastal Management,
#' 60, 19-30. doi: 10.1016/j.ocecoaman.2011.12.016
#' 
#' The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces of the World. GIS layers
#' developed by The Nature Conservancy with multiple partners, combined from Spalding et al. (2007)
#' and Spalding et al. (2012). Cambridge (UK): The Nature Conservancy.
#' 
#' Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice Kottelat, Nina Bogutskaya,
#' Brian Coad, Nick Mandrak, Salvador Contreras Balderas, William Bussing, Melanie L. J. Stiassny,
#' Paul Skelton, Gerald R. Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf,
#' James Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel, Eric Wikramanayake,
#' David Olson, Hugo L. López, Roberto E. Reis, John G. Lundberg, Mark H. Sabaj Pérez,
#' Paulo Petry, Freshwater Ecoregions of the World: A New Map of Biogeographic Units for
#' Freshwater Biodiversity Conservation, BioScience, Volume 58, Issue 5, May 2008,
#' Pages 403–414. doi: 10.1641/B580507
#' 
#' Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Terra - CRAN
#' @seealso
#' For more information on the original code and methods, check Hagen, Oskar et al. (2019), Data
#' from: Mountain building, climate cooling and the richness of cold-adapted plants in the northern
#' hemisphere, Dryad, Dataset, https://doi.org/10.5061/dryad.0ff6b04.
#' @example inst/examples/get_range_help.R
#' @importFrom rnaturalearth ne_countries
#' @importFrom methods is
#' @importFrom terra vect crds intersect simplifyGeom buffer rast disagg aggregate rasterize crop
#' @importFrom FNN knn.dist
#' @importFrom stats kmeans
#' @importFrom mclust Mclust mclustBIC
#' @importFrom ClusterR KMeans_rcpp
#' @export
get_range <- function (occ_coord = NULL, 
                       bioreg = NULL, 
                       bioreg_name = NULL, 
                       degrees_outlier = 5,
                       clustered_points_outlier = 3,
                       buffer_width_point = 4, 
                       buffer_increment_point_line = 0.5, 
                       buffer_width_polygon = 4,
                       dir_temp = tempdir(),
                       raster = TRUE,
                       res = 0.01,
                       verbose = TRUE){

  ### =========================================================================
  ### Object conditions + remove duplicates
  ### =========================================================================

  # Bioreg and convert to sf
  if (!class(bioreg)[1] %in% c("SpatialPolygonsDataFrame", "SpatVector", "sf")) {
     stop("Wrong 'bioreg' class (not a spatial object)...")
  }
  if (class(bioreg)[1] %in% c("SpatialPolygonsDataFrame", "sf")) {
    bioreg <- terra::vect(bioreg)
  }

  # Bioreg name
  if (names(bioreg)[1] %in% "CLARA") {
    bioreg_name <- "EcoRegion"
  } 
  if (methods::is(bioreg_name, "NULL")) {
    stop("Name of the desired ecoregion level/category ('bioreg_name' parameter) is missing, please provide one...")
  } 

  # occ_coord
  if (!methods::is(occ_coord, "data.frame")) {
    stop("'occ_coord' is not a data.frame...")
  } 
  if (!any(names(occ_coord) %in% "decimalLongitude")) {
    stop("Longitute/Latitude columns wrongly defined...")
  }
  
  # get sp.name, and stop if not unique
  sp.name <- unique(occ_coord$input_search)
  if (length(sp.name) > 1) {
    stop("More than one species in the input data.frame... \n
         gbif.range is meant to be used for one species at a time. \n
         Only consider multiple species (with a common ID) if they are closely related and share similar ecological characteristics.")
  }

  # Remove duplicates
  w.col <- c("decimalLongitude","decimalLatitude")
  occ.coord.k <- occ_coord
  occ_coord[, w.col] <- round(occ_coord[,w.col], 4)
  occ_coord <- occ_coord[!duplicated(occ_coord[, w.col]), ]
  occ_coord <- terra::vect(occ_coord, geom = c("decimalLongitude","decimalLatitude"), crs = "epsg:4326")

  ### =========================================================================
  ### Check if sufficient data
  ### =========================================================================
  
  # Check if there sufficient species & if not, make an entry in the log-file and end the function
  if (nrow(occ_coord) <= clustered_points_outlier +1 ){
    stop("Too few occurences!")
  } 
  
  if (verbose){
    cat("## Start of computation for species: ", sp.name, " ###", "\n") 
  }
  
  ### =========================================================================
  ### Identify outliers
  ### =========================================================================
  
  # Create distance matrix...
  mat.dist <- as.matrix(FNN::knn.dist(terra::crds(occ_coord), k = clustered_points_outlier))
  
  # Mark outliers
  cond <- apply(mat.dist, 1, function(x) x[clustered_points_outlier]) > degrees_outlier
  
  # Print info
  if (verbose){
    cat(paste0(sum(cond), " outlier's from " ,nrow(occ_coord), " | proportion from total points: ",
      round((sum(cond) / nrow(occ_coord)) * 100,0), "%\n"))
  }
  
  # Remove outliers
  occ.coord.mod <- occ_coord[!cond, ]

  # Stop if too many outliers in data set
  if(nrow(occ.coord.mod) == 0){
    stop('Too few occurrences within outlier threshold!')
  } 

  ### =========================================================================
  ### Define those polygons
  ### =========================================================================
  
  # Set number of ecoregions
  ovo.coord.mod <- terra::intersect(bioreg,occ.coord.mod)
  uniq <- levels(factor(ovo.coord.mod[[bioreg_name]][[1]]))

  # Handling NA error
  if (length(ovo.coord.mod) == 0) {
    stop("No overlap of ecoregion info, please use another 'bioreg_name'")
  }
  
  # Loop over bioregions
  SP.dist <- list()
  for(g in 1:length(uniq)) {
    
    # Print
    if (verbose){
      cat('bioregion', g, ' of ',length(uniq),": ",uniq[g], '\n')
    }

    # NAs or not
    q1 <- as.data.frame(bioreg)[, bioreg_name] == uniq[g]
    q1[is.na(q1)] <- FALSE

    # Continue
    tmp <- terra::simplifyGeom(bioreg[q1, ], tolerance = 0.001, preserveTopology = TRUE)
    a <- ovo.coord.mod[ovo.coord.mod[[bioreg_name]][[1]] == uniq[g]]
    
    if (length(a) < 3) {
      k <- 1
      cluster.k <- stats::kmeans(terra::crds(a), k)
      cluster.k$clusters <- cluster.k$cluster 

    } else {

      if (all(terra::crds(a)[, 1] == terra::crds(a)[, 2])) {
        terra::crds(a)[, 2] <- terra::crds(a)[, 2] + 0.00001
      }
      
      # Determine number of clusters
      m.clust <- mclust::Mclust(terra::crds(a) + 1000, verbose = FALSE)
      
      # k = number of clusters
      k <- m.clust$G 
      
      # Reduce k if necessary so that KMeans_rcpp() will run
      while (k > length(a)-2) {k <- k-1} 
      if (k == 0) {k <- 1}
      
      cluster.k <- ClusterR::KMeans_rcpp(terra::crds(a), k, num_init = 20, initializer = 'random')
      
      while (length(unique(cluster.k$clusters)) < k) {
        k <- k-1
        cluster.k <- ClusterR::KMeans_rcpp(terra::crds(a), k, num_init = 20, initializer = 'random')
      }
      
    }
   
    polygons.list <- list() 
    for (i in 1:k)
    {
      # kmeans (with number of clusters from mcluster)
      a.temp <- a[cluster.k$clusters == i]

      # Generate polygon
      my.shpe <- conv_function(sp_coord = a.temp,
                            bwp = buffer_width_point,
                            bipl = buffer_increment_point_line,
                            bwpo = buffer_width_polygon,
                            temp_dir = dir_temp,
                            g = g)
      
      # Intersect polygon with ecoregion (zero buffer to avoid error)
      b1 <- terra::buffer(my.shpe, width = 0)
      b2 <- terra::buffer(tmp, width = 0)
      polygons.list[[i]] <- terra::intersect(b1, b2)
    }  
    
    SP.dist[[g]] <- do.call("rbind", polygons.list)
  } 
  
  lala <- SP.dist[!is.na(SP.dist)]
  
  ### =========================================================================
  ### Check and return output
  ### =========================================================================
  
  if (!dir.exists(dir_temp)) {
    unlink(dir_temp, recursive = TRUE)
  }
  
  if (length(lala) == 0) {
    stop('No occurrences within Bioregions. Empty raster produced...')
  }
  shp.species <- do.call("rbind", lala)

  # Convert to raster or not
  if (raster) {
    res.use = 1 / res
    ras.res <- terra::rast(terra::disagg(terra::rast(), res.use))
    sp.range.u <- terra::aggregate(shp.species)
    ras <- terra::rasterize(sp.range.u, ras.res)
    shp.species <- terra::crop(ras, sp.range.u)
  }
  
  # Final print
  if (verbose){
    cat("## End of computation for species: ", sp.name, " ###", "\n")
  }

  # Out
  result <- list(init.args = list(
                  occ_coord = occ.coord.k,
                  bioreg = bioreg,
                  bioreg_name = bioreg_name, 
                  degrees_outlier = degrees_outlier,
                  clustered_points_outlier = clustered_points_outlier,
                  buffer_width_point = buffer_width_point,
                  buffer_increment_point_line = buffer_increment_point_line,
                  buffer_width_polygon = buffer_width_polygon,
                  dir_temp = dir_temp,
                  raster = TRUE,
                  res = res
                ),
                rangeOutput = shp.species)

  return(result)
}
