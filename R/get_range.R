### ==================================================================
### get_range
### ==================================================================
#' Create a species range map based on \code{get_gbif()} object
#'
#' This function estimates species ranges from occurrence data (GBIF or
#' else) and ecoregions [see \code{make_ecoregion()} or \code{bioreg_list}].
#' 
#' @param occ_coord Object of class \code{getGBIF} [see \code{get_gbif()}]
#' or a \code{data.frame} containing two columns named 'decimalLongitude'
#' and 'decimalLatitude'.
#' @param bioreg  Object of class \code{SpatialPolygonsDataFrame},
#' \code{SpatVector} or \code{sf} containing different ecoregions (WGS84).
#' Define the range extent and ecoregions. Note that this
#' parameter may be fed with an external, generated
#' [see \code{make_ecoregion()}] or downloaded ecoregion shapefile
#' [see \code{read_bioreg()}]. Four ecoregions can be downloaded
#' with the library: 'eco_terra', 'eco_marine', 'eco_hd_marine'
#' and 'eco_fresh' (see \code{bioreg_list} and details below).
#' @param bioreg_name Character. One ecoregion level/category name from the
#' \code{bioreg} parameter polygon must be supplied, e.g., very detailed level
#' of eco_terra' is "ECO_NAME". Note that default applies if a
#' \code{make_ecoregion()} polygon is provided as a \code{bioreg} parameter.
#' @param degrees_outlier Numeric. Distance threshold (degrees). Points whose 
#' `\code{clust_pts_outlier}-th nearest neighbour exceeds this distance are
#' classified as outliers (default: 5°).
#' @param clust_pts_outlier Numeric. k-NN order for outlier detection. Points 
#' must have ≥\code{k} nearest neighbors within \code{degrees_outlier} distance
#' to be retained (default: 4).
#' @param buffer_width_point Numeric. Buffer (in degrees) which will be applied
#' around single observations.
#' @param buff_incrmt_pt_line Numeric. How much should the buffer be increased
#' for each point on a line.
#' @param buffer_width_polygon Numeric. Buffer (in degrees) which will be
#' applied around distribution polygons (for each ecoregion).
#' @param dir_temp Character. Where should the temporary text file for the
#' convex hull be saved? (text file will be deleted again). Default value
#' is \code{tempdir()}.
#' @param raster Logical. Should the output be a unified raster?
#' Default is \code{TRUE}
#' @param res Numeric. If \code{raster = TRUE}, resolution of the output in
#' degrees (1° = ~111 km at the equator). Default is 0.1 (~11.1 km). It is
#' important to note that the highest achievable resolution of the output will
#' depend on its \code{bioreg} precision, e.g., a species range output can reach
#' the same resolution of the rasters used to create a \code{make_ecoregion()}
#' object.
#' @param verbose Logical. Should progession be printed?
#' @details
#' The function implements a four-step species range mapping process: 
#' 
#' **Step 1 - Outlier filtering and ecoregion assignement**: Outliers
#' are removed from occurrence records using k-nearest neighbor (k-NN)
#' distances. Points whose distance to their \code{clust_pts_outlier}-th
#' nearest neighbor exceeds the  \code{degrees_outlier} threshold (default:
#' 5°) are excluded, retaining only well-supported clusters (default: ≥4
#' points within 5°). Then, non-outlier points are spatially intersected 
#' with ecoregions (\code{bioreg}, specified via \code{bioreg_name}) to
#' identify occupied bioregions.
#' 
#' **Step 2 - Clustering**: Within each occupied ecoregion, points are
#' clustered using Gaussian mixture modeling (\code{mclust::Mclust}) to
#' determine the optimal number of clusters, followed by k-means clustering 
#' (\code{ClusterR::KMeans_rcpp}). For <3 points, a single cluster is assigned 
#' via basic k-means; collinear points receive minimal jittering.
#' 
#' **Step 3 - Convex hull**: For each cluster within an ecoregion, 
#' \code{conv_function()} generates a buffered convex hull. Special cases
#' include: (1) single points receive circular buffers
#' (\code{buffer_width_point}, default: 4°); (2) collinear points
#' (suggesting transects) receive incrementally widening buffers along
#' the line (\code{buff_incrmt_pt_line}, default: 0.5° per additional
#' point); and (3) standard clusters receive polygon expansion
#' (\code{buffer_width_polygon}, default: 4°).
#'
#' **Step 4 - Ecological intersection**: Each cluster-derived polygon is 
#' intersected with its parent ecoregion boundary (after zero-width buffering 
#' to ensure topological validity). Per-ecoregion cluster polygons are merged, 
#' then combined across all occupied ecoregions to produce the final species 
#' range (\code{SpatVector} or \code{SpatRaster} at \code{res = 0.1°} if 
#' \code{raster = TRUE}).
#''
#' If there are too many records, \code{get_range()} can process a sub-sample
#' of species observations to speed up polygon creation or avoid potential RAM
#' issues.
#' 
#' Ecoregions represent large, geographically distinct areas containing 
#' characteristic assemblages of species, dynamics, and environmental conditions 
#' (https://en.wikipedia.org/wiki/Ecoregion). Download-ready ecoregion datasets 
#' include:
#'
#' (1) 'eco_terra' is for terrestrial species. Has three different levels:
#' "ECO_NAME", "WWF_MHTNAM" and "WWF_REALM2".
#'
#' (2) 'eco_fresh' is for freshwater species. Has only one: "ECOREGION".
#' 
#' (3) 'eco_marine' and 'eco_hd_marine' is for marine species. It contains
#' three distinct levels: "ECOREGION", "PROVINCE" and "REALM".
#'
#' For marine species, 'eco_terra' may also be used if the user wants to
#' represent the terrestrial range of species that also partially settle
#' on mainland. For fresh water species, same may be done if the user
#' considers that terrestrial ecoregions should be more representative of
#' the species ecology.
#' @return An object of class \code{getRange} with two fields:
#' \code{init.args} (parameters and data employed) and
#' \code{rangeOutput} (object of class \code{SpatVector} or
#' \code{SpatRaster} depending on what the user set as \code{raster}
#' parameter).
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
#' of life on Earth. Bioscience 51(11):933-938. doi: 10.1641/0006-3568(2001)051
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
#' @seealso
#' For more information on the original code and methods, check Hagen et al.
#' (2019), Data from: Mountain building, climate cooling and the richness
#' of cold-adapted plants in the northern hemisphere, Dryad, Dataset,
#' https://doi.org/10.5061/dryad.0ff6b04.
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
                       bioreg = NULL, 
                       bioreg_name = NULL, 
                       degrees_outlier = 5,
                       clust_pts_outlier = 4,
                       buffer_width_point = 4, 
                       buff_incrmt_pt_line = 0.5, 
                       buffer_width_polygon = 4,
                       dir_temp = tempdir(),
                       raster = TRUE,
                       res = 0.1,
                       verbose = TRUE){

  ######################################################
  ### Stop messages
  ######################################################


  check_numeric(degrees_outlier, "degrees_outlier")
  check_numeric(clust_pts_outlier, "clust_pts_outlier")
  check_numeric(buffer_width_point, "buffer_width_point")
  check_numeric(buff_incrmt_pt_line, "buff_incrmt_pt_line")
  check_numeric(buffer_width_polygon, "buffer_width_polygon")
  check_logical(raster, "raster")
  check_numeric(res, "res")
  check_logical(verbose, "verbose")

  # Spatial class
  spatial.class <- c("SpatialPolygonsDataFrame",
    "SpatialPolygons","sf", "SpatVector")

  # Bioreg and convert to sf
  if (!class(bioreg)[1] %in% spatial.class) {
     stop("Wrong 'bioreg' class (not a spatial object)...")
  }
  if (class(bioreg)[1] %in% spatial.class[c(1,2)]) {
    bioreg <- terra::vect(bioreg)
  }

  # Bioreg name
  if (names(bioreg)[1] %in% "CLARA") {
    bioreg_name <- "EcoRegion"
  } 
  if (methods::is(bioreg_name, "NULL")) {
    stop("Name of the desired ecoregion level/category
      ('bioreg_name' parameter) is missing, please provide one...")
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
    stop("Too few occurences!")
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
    x = bioreg,
    y = occ.coord.mod
  )
  uniq <- levels(factor(ovo.coord.mod[[bioreg_name]][[1]]))

  # Handling NA error
  if (length(ovo.coord.mod) == 0) {
    stop("No overlap of ecoregion info, please use another 'bioreg_name'")
  }
  
  # Loop over bioregions
  SP.dist <- list()
  for (g in seq_along(uniq)) {
    
    # Print
    if (verbose){
      cat('bioregion', g, ' of ',length(uniq),": ",uniq[g], '\n')
    }

    # NAs or not
    q1 <- as.data.frame(bioreg)[, bioreg_name] == uniq[g]
    q1[is.na(q1)] <- FALSE

    # Continue
    tmp <- terra::simplifyGeom(
      x = bioreg[q1, ],
      tolerance = 0.001,
      preserveTopology = TRUE
    )
    a <- ovo.coord.mod[ovo.coord.mod[[bioreg_name]][[1]] == uniq[g]]
    
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
        bwp = buffer_width_point,
        bipl = buff_incrmt_pt_line,
        bwpo = buffer_width_polygon,
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
    stop('No occurrences within Bioregions. Empty raster produced...')
  }
  shp.species <- do.call("rbind", lala)

  # Convert to raster or not
  if (raster) {
    res.use <- 1 / res
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
  result <- getRange$new()
  result$init.args <- list(occ_coord = occ.coord.k,
                           bioreg = bioreg,
                           bioreg_name = bioreg_name, 
                           degrees_outlier = degrees_outlier,
                           clust_pts_outlier = clust_pts_outlier,
                           buffer_width_point = buffer_width_point,
                           buff_incrmt_pt_line = buff_incrmt_pt_line,
                           buffer_width_polygon = buffer_width_polygon,
                           dir_temp = dir_temp,
                           raster = TRUE,
                           res = res)
  result$rangeOutput <- shp.species

  return(result)
}
