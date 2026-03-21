### ==================================================================
### get_range
### ==================================================================
#' Create Species Range Maps from Occurrences and Ecoregions
#'
#' Estimate ecologically informed species ranges from occurrence data and an
#' ecoregion layer. The workflow combines outlier filtering, clustering, convex
#' hull construction, and intersection with occupied ecoregions.
#' 
#' @param occ_coord A \code{getGBIF} object returned by \code{get_gbif()} or a
#' \code{data.frame} containing the columns \code{decimalLongitude} and
#' \code{decimalLatitude}.
#' @param ecoreg Spatial ecoregion layer in WGS84. Accepted classes are
#' \code{SpatialPolygonsDataFrame}, \code{SpatVector}, and \code{sf}. This can
#' be a downloaded layer from \code{read_ecoreg()} or a custom layer created by
#' \code{make_ecoreg()}.
#' @param ecoreg_name Character string naming the field in \code{ecoreg} that
#' defines ecoregion categories. For \code{eco_terra}, for example,
#' \code{"ECO_NAME"} is the most detailed level. When \code{ecoreg} comes from
#' \code{make_ecoreg()}, \code{"EcoRegion"} is used automatically.
#' @param degrees_outlier Numeric distance threshold in degrees. Points whose
#' \code{clust_pts_outlier}-th nearest neighbour lies beyond this distance are
#' treated as outliers. Default is \code{5}.
#' @param clust_pts_outlier Numeric k-nearest-neighbour order used for outlier
#' detection. Default is \code{4}.
#' @param buff_width_point Numeric buffer width in degrees for isolated points.
#' @param buff_incrmt_pts_line Numeric increment in buffer width for linear
#' clusters.
#' @param buff_width_polygon Numeric buffer width in degrees applied to convex
#' hull polygons.
#' @param dir_temp Character string giving the directory used for temporary
#' convex-hull files. Defaults to \code{tempdir()}.
#' @param raster Logical. Should the final output be rasterized? Default is
#' \code{TRUE}.
#' @param format Output format used when \code{raster = FALSE}. Choose between
#' \code{"SpatVector"} (default) and \code{"sf"}.
#' @param res Numeric output resolution in degrees when \code{raster = TRUE}.
#' Default is \code{0.1} (about 11.1 km at the equator). The achievable
#' resolution is constrained by the spatial precision of \code{ecoreg}.
#' @param verbose Logical. Should progress messages be printed?
#' @details The function follows four main steps.
#'
#' First, occurrence points are filtered for spatial outliers using
#' k-nearest-neighbour distances and are assigned to ecoregions.
#'
#' Second, points within occupied ecoregions are grouped into clusters using a
#' combination of Gaussian mixture modeling and k-means clustering.
#'
#' Third, each cluster is converted into a polygon. Single points receive
#' circular buffers, collinear clusters receive line-based buffers, and other
#' clusters receive buffered convex hulls through \code{conv_function()}.
#'
#' Fourth, cluster polygons are intersected with their parent ecoregions and
#' merged into a final range layer.
#'
#' Download-ready ecoregion datasets include \code{eco_terra} for terrestrial
#' species, \code{eco_fresh} for freshwater species, and \code{eco_marine} or
#' \code{eco_hd_marine} for marine species.
#' @return An object of class \code{getRange} with two fields:
#' \code{init.args}, containing the arguments and data used to build the map,
#' and \code{rangeOutput}, containing the resulting \code{SpatVector},
#' \code{sf}, or \code{SpatRaster} object.
#' @references
#' Oskar Hagen, Lisa Vaterlaus, Camille Albouy, Andrew Brown, Flurin Leugger,
#' Renske E. Onstein, Charles Novaes de Santana, Christopher R. Scotese,
#' Loïc Pellissier. (2019) Mountain building, climate cooling and the
#' richness of cold-adapted plants in the Northern Hemisphere. Journal of
#' Biogeography. doi: 10.1111/jbi.13653
#' 
#' Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J. S.,
#' ... &  Pellissier, L. (2022). An integrated high resolution mapping shows
#' congruent biodiversity patterns of Fagales and Pinales. New Phytologist,
#' 235(2), 759-772 10.1111/nph.18158
#' 
#' Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
#' Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I.,
#' Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F.,
#' Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao,
#' P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map
#' of life on Earth. BioScience 51(11):933-938. doi: 10.1641/0006-3568(2001)051
#' 
#' The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
#' Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
#' Units. GIS layers developed by The Nature Conservancy with multiple partners,
#' combined from Olson et al. (2001), Bailey 1995 and Wiken 1986. Cambridge
#' (UK): The Nature Conservancy.
#' 
#' Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A.
#' Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana,
#' Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer Molnar, Cheri
#' A. Recchia, James Robertson, Marine Ecoregions of the World: A
#' Bioregionalization of Coastal and Shelf Areas, BioScience, Volume 57,
#' Issue 7, July 2007, Pages 573–583. doi: 10.1641/B570707
#' 
#' Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
#' Pelagic provinces of the world: a biogeographic classification of the
#' world’s surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
#' doi: 10.1016/j.ocecoaman.2011.12.016
#' 
#' The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
#' of the World. GIS layers developed by The Nature Conservancy with multiple
#' partners, combined from Spalding et al. (2007) and Spalding et al. (2012).
#' Cambridge (UK): The Nature Conservancy.
#' 
#' Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice
#' Kottelat, Nina Bogutskaya, Brian Coad, Nick Mandrak, Salvador Contreras
#' Balderas, William Bussing, Melanie L. J. Stiassny, Paul Skelton, Gerald R.
#' Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf, James
#' Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel, Eric
#' Wikramanayake, David Olson, Hugo L. López, Roberto E. Reis, John G.
#' Lundberg, Mark H. Sabaj Pérez, Paulo Petry, Freshwater Ecoregions of
#' the World: A New Map of Biogeographic Units for Freshwater Biodiversity
#' Conservation, BioScience, Volume 58, Issue 5, May 2008, Pages 403–414.
#' doi: 10.1641/B580507
#' 
#' Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7."
#' (2022). Terra - CRAN
#' @seealso \code{read_ecoreg()}, \code{make_ecoreg()}, and \code{cv_range()}.
#' @example inst/examples/get_range_help.R
#' @importFrom rnaturalearth ne_countries
#' @importFrom methods is new
#' @importFrom terra vect crds intersect simplifyGeom buffer
#' rast disagg aggregate rasterize crop
#' @importFrom FNN knn.dist
#' @importFrom stats kmeans
#' @importFrom mclust Mclust mclustBIC
#' @importFrom ClusterR KMeans_rcpp
#' @export
get_range <- function (occ_coord = NULL, 
                       ecoreg = NULL, 
                       ecoreg_name = NULL, 
                       degrees_outlier = 5,
                       clust_pts_outlier = 4,
                       buff_width_point = 4, 
                       buff_incrmt_pts_line = 0.5, 
                       buff_width_polygon = 4,
                       dir_temp = tempdir(),
                       raster = TRUE,
                       format = "SpatVector",
                       res = 0.1,
                       verbose = TRUE){

  ######################################################
  ### Stop messages
  ######################################################


  check_numeric(degrees_outlier, "degrees_outlier")
  check_numeric(clust_pts_outlier, "clust_pts_outlier")
  check_numeric(buff_width_point, "buff_width_point")
  check_numeric(buff_incrmt_pts_line, "buff_incrmt_pts_line")
  check_numeric(buff_width_polygon, "buff_width_polygon")
  check_logical(raster, "raster")
  check_numeric(res, "res")
  check_logical(verbose, "verbose")

  # Spatial class
  spatial.class <- c("SpatialPolygonsDataFrame",
    "SpatialPolygons","sf", "SpatVector")

  # Ecoreg and convert to sf
  if (!class(ecoreg)[1] %in% spatial.class) {
     stop("Wrong 'ecoreg' class (not a spatial object)...")
  }
  if (class(ecoreg)[1] %in% spatial.class[1:3]) {
    ecoreg <- terra::vect(ecoreg)
  }

  # Ecoreg name
  if (names(ecoreg)[1] %in% "CLARA") {
    ecoreg_name <- "EcoRegion"
  } 
  if (methods::is(ecoreg_name, "NULL")) {
    stop(
      paste(
        "Name of the desired ecoregion level/category,
        ('ecoreg_name' parameter)",
        "is missing, please provide one..."
      )
)
  } 

  # occ_coord
  if (!methods::is(occ_coord, "data.frame")) {
    stop("'occ_coord' is not a data.frame...")
  } 
  if (!any(names(occ_coord) %in% "decimalLongitude")) {
    stop("Longitude/Latitude columns wrongly defined...")
  }
  
  # get sp.name, and stop if not unique
  sp.name <- unique(occ_coord$input_search)
  if (length(sp.name) > 1) {
    stop(
      paste(
      "More than one species in the input data.frame...",
      "gbif.range is meant to be used for one species at a time.",
      "Only consider multiple species (with a common ID) if they are",
      "closely related and share similar ecological characteristics."
    ))
  }

  # Remove duplicates
  w.col <- c("decimalLongitude","decimalLatitude")
  occ.coord.k <- occ_coord
  occ_coord[, w.col] <- round(occ_coord[,w.col], 4)
  occ_coord <- occ_coord[!duplicated(occ_coord[, w.col]), ]
  occ_coord <- terra::vect(
    x = as.data.frame(occ_coord),
    geom = c("decimalLongitude","decimalLatitude"),
    crs = "epsg:4326"
  )


  ### =========================================================================
  ### Check if sufficient data
  ### =========================================================================

  
  # Check if there sufficient species & if not, make an entry in the
  # log-file and end the function
  if (nrow(occ_coord) <= clust_pts_outlier +1){
    stop("Too few occurrences!")
  } 
  
  if (verbose){
    cat("## Start of computation for species: ", sp.name, " ###", "\n") 
  }
  

  ### =========================================================================
  ### Identify outliers
  ### =========================================================================
  

  # Create distance matrix...
  mat.dist <- as.matrix(
    FNN::knn.dist(
      data = terra::crds(occ_coord),
      k = clust_pts_outlier
    )
  )
  
  # Mark outliers
  cond <- apply(mat.dist, 1,
    function(x) x[clust_pts_outlier]) > degrees_outlier
  
  # Print info
  if (verbose){
    cat(
      paste0(sum(cond)," outlier's from " ,
        nrow(occ_coord), " | proportion from total points: ",
        round((sum(cond) / nrow(occ_coord)) * 100,0), "%\n"
      )
    )
  }
  
  # Remove outliers
  occ.coord.mod <- occ_coord[!cond, ]

  # Stop if too many outliers in data set
  if (nrow(occ.coord.mod) == 0){
    stop('Too few occurrences within outlier threshold!')
  } 


  ### =========================================================================
  ### Define those polygons
  ### =========================================================================
  

  # Set number of ecoregions
  ovo.coord.mod <- terra::intersect(
    x = ecoreg,
    y = occ.coord.mod
  )
  uniq <- levels(factor(ovo.coord.mod[[ecoreg_name]][[1]]))

  # Handling NA error
  if (length(ovo.coord.mod) == 0) {
    stop("No overlap of ecoregion info, please use another 'ecoreg_name'")
  }
  
  # Loop over ecoregions
  SP.dist <- list()
  for (g in seq_along(uniq)) {
    
    # Print
    if (verbose){
      cat('ecoregion', g, ' of ',length(uniq),": ",uniq[g], '\n')
    }

    # Handle NAsy
    q1 <- as.data.frame(ecoreg)[, ecoreg_name] == uniq[g]
    q1[is.na(q1)] <- FALSE

    # Continue
    tmp <- terra::simplifyGeom(
      x = ecoreg[q1, ],
      tolerance = 0.001,
      preserveTopology = TRUE
    )
    a <- ovo.coord.mod[ovo.coord.mod[[ecoreg_name]][[1]] == uniq[g]]
    
    if (length(a) < 3) {
      k <- 1
      cluster.k <- stats::kmeans(
        x = terra::crds(a),
        centers = k
      )
      cluster.k$clusters <- cluster.k$cluster 

    } else {

      if (all(terra::crds(a)[, 1] == terra::crds(a)[, 2])) {
        terra::crds(a)[, 2] <- terra::crds(a)[, 2] + 0.00001
      }
      
      # Determine number of clusters
      m.clust <- mclust::Mclust(
        data = terra::crds(a) + 1000,
        verbose = FALSE
      )
      
      # k = number of clusters
      k <- m.clust$G 
      
      # Reduce k if necessary so that KMeans_rcpp() will run
      while (k > length(a)-2) {k <- k-1} 
      if (k == 0) {k <- 1}
      
      cluster.k <- ClusterR::KMeans_rcpp(
        data = terra::crds(a),
        clusters = k,
        num_init = 20,
        initializer = 'random'
      )
      
      while (length(unique(cluster.k$clusters)) < k) {
        k <- k-1
        cluster.k <- ClusterR::KMeans_rcpp(
          data = terra::crds(a),
          clusters = k,
          num_init = 20,
          initializer = 'random'
        )
      }
    }
   
    polygons.list <- list() 
    for (i in seq_len(k))
    {
      # kmeans (with number of clusters from mcluster)
      a.temp <- a[cluster.k$clusters == i]

      # Generate polygon
      my.shpe <- conv_function(
        sp_coord = a.temp,
        bwp = buff_width_point,
        bipl = buff_incrmt_pts_line,
        bwpo = buff_width_polygon,
        temp_dir = dir_temp,
        g = g
      )
      
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
    stop('No occurrences within ecoregions. Empty raster produced...')
  }
  shp.species <- do.call("rbind", lala)

  # Convert to raster or not
  if (raster) {
    res.use <- 1 / res
    ras.res <- terra::rast(terra::disagg(terra::rast(), res.use))
    sp.range.u <- terra::aggregate(shp.species)
    ras <- terra::rasterize(sp.range.u, ras.res)
    shp.species <- terra::crop(ras, sp.range.u)
    names(shp.species) <- occ_coord$input_search[1]

  } else if (format == "sf") {
    shp.species <- sf::st_as_sf(shp.species)
  }
  
  # Final print
  if (verbose){
    cat("## End of computation for species: ", sp.name, " ###", "\n")
  }

  # Out
  result <- getRange$new()
  result$init.args <- list(occ_coord = occ.coord.k,
                           ecoreg = ecoreg,
                           ecoreg_name = ecoreg_name, 
                           degrees_outlier = degrees_outlier,
                           clust_pts_outlier = clust_pts_outlier,
                           buff_width_point = buff_width_point,
                           buff_incrmt_pts_line = buff_incrmt_pts_line,
                           buff_width_polygon = buff_width_polygon,
                           dir_temp = dir_temp,
                           raster = TRUE,
                           res = res)
  result$rangeOutput <- shp.species

  return(result)
}
