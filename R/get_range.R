### ==================================================================
### get_range
### ==================================================================
#' Create Species Range Maps from Occurrences and Ecoregions
#'
#' Estimate ecologically informed species ranges from occurrence data and an
#' ecoregion layer. The workflow combines outlier filtering, clustering, convex
#' hull construction, and intersection with occupied ecoregions.
#' 
#' @param occ_coord \code{getGBIF} object returned by \code{get_gbif()} or a
#' \code{data.frame} containing the columns \code{decimalLongitude} and
#' \code{decimalLatitude}.
#' @param ecoreg Spatial ecoregion layer in WGS84. Accepted classes are
#' \code{SpatialPolygonsDataFrame}, \code{SpatVector}, and \code{sf}. This can
#' be a downloaded layer from \code{read_ecoreg()} or a custom layer created by
#' \code{make_ecoreg()}.
#' @param ecoreg_name Character. String naming the field in \code{ecoreg} that
#' defines ecoregion categories. For \code{eco_terra}, for example,
#' \code{"ECO_NAME"} is the most detailed level. When \code{ecoreg} comes from
#' \code{make_ecoreg()}, \code{"EcoRegion"} is used automatically.
#' @param degrees_outlier Numeric. Distance threshold in degrees. Points whose
#' \code{clust_pts_outlier}-th nearest neighbour lies beyond this distance are
#' treated as outliers. Default is \code{5}.
#' @param clust_pts_outlier Numeric. k-nearest-neighbour order used for outlier
#' detection. Default is \code{4}.
#' @param buff_width_point Numeric. Buffer width in degrees for isolated points.
#' @param buff_incrmt_pts_line Numeric. Increment in buffer width for linear
#' clusters.
#' @param buff_width_polygon Numeric. Buffer width in degrees applied to convex
#' hull polygons.
#' @param dir_temp Character. String giving the parent directory used for
#' temporary convex-hull files. A uniquely named sub-directory is created
#' inside it and removed when the function exits; \code{dir_temp} itself is
#' never deleted. Defaults to \code{tempdir()}.
#' @param format Character. Output format for the range map. One of
#' \code{"SpatVector"} (default), \code{"sf"}, or \code{"SpatRaster"}.
#' \code{"SpatRaster"} rasterizes the range at the resolution set by
#' \code{res}.
#' @param res Numeric. Output resolution in degrees when \code{format = "SpatRaster"}.
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
#'
#' For \code{format = "SpatVector"} and \code{format = "sf"}, the range is
#' returned as \strong{one feature per occupied ecoregion} rather than a
#' single merged polygon. Features may be multipart where a species occupies
#' several disjoint patches of the same ecoregion; \code{terra::disagg()}
#' separates them. Because ecoregions do not overlap, the features do not
#' overlap either, so their areas are additive and
#' \code{sum(terra::expanse(x))} gives the total range size. Each feature
#' carries the \code{ecoreg_name} field it was clipped to plus a
#' \code{species} column. Use \code{merge_range()} to dissolve the features
#' into one polygon, or \code{format = "SpatRaster"} for a gridded range.
#' @references
#' Hagen, O., Vaterlaus, L., Albouy, C., Brown, A., Leugger, F., Onstein,
#' R. E., Novaes de Santana, C., Scotese, C. R., & Pellissier, L. (2019).
#' Mountain building, climate cooling and the richness of cold-adapted
#' plants in the Northern Hemisphere. Journal of Biogeography.
#' \doi{10.1111/jbi.13653}
#'
#' Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J. S.,
#' ... & Pellissier, L. (2022). An integrated high resolution mapping shows
#' congruent biodiversity patterns of Fagales and Pinales. New Phytologist,
#' 235(2), 759-772. \doi{10.1111/nph.18158}
#'
#' Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
#' Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I.,
#' Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F.,
#' Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao,
#' P., Kassem, K. R. (2001). Terrestrial ecoregions of the world: a new map
#' of life on Earth. BioScience, 51(11), 933-938.
#' \doi{10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2}
#'
#' The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
#' Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
#' Units. GIS layers developed by The Nature Conservancy with multiple
#' partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986.
#' Cambridge (UK): The Nature Conservancy.
#'
#' Spalding, M. D., Fox, H. E., Allen, G. R., Davidson, N., Ferda\enc{ñ}{n}a, Z. A.,
#' Finlayson, M., Halpern, B. S., Jorge, M. A., Lombana, A., Lourie, S. A.,
#' Martin, K. D., McManus, E., Molnar, J., Recchia, C. A., Robertson, J.
#' (2007). Marine Ecoregions of the World: A Bioregionalization of Coastal
#' and Shelf Areas. BioScience, 57(7), 573-583. \doi{10.1641/B570707}
#'
#' Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
#' Pelagic provinces of the world: a biogeographic classification of the
#' world's surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
#' \doi{10.1016/j.ocecoaman.2011.12.016}
#'
#' The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
#' of the World. GIS layers developed by The Nature Conservancy with multiple
#' partners, combined from Spalding et al. (2007) and Spalding et al. (2012).
#' Cambridge (UK): The Nature Conservancy.
#'
#' Abell, R., Thieme, M. L., Revenga, C., Bryer, M., Kottelat, M.,
#' Bogutskaya, N., Coad, B., Mandrak, N., Contreras Balderas, S., Bussing,
#' W., Stiassny, M. L. J., Skelton, P., Allen, G. R., Unmack, P., Naseka,
#' A., Ng, R., Sindorf, N., Robertson, J., Armijo, E., Higgins, J. V.,
#' Heibel, T. J., Wikramanayake, E., Olson, D., Lopez, H. L., Reis, R. E.,
#' Lundberg, J. G., Sabaj Perez, M. H., Petry, P. (2008). Freshwater
#' Ecoregions of the World: A New Map of Biogeographic Units for Freshwater
#' Biodiversity Conservation. BioScience, 58(5), 403-414.
#' \doi{10.1641/B580507}
#'
#' Hijmans, R. J. (2022). terra: Spatial Data Analysis. R package version
#' 1.6-7. \url{https://CRAN.R-project.org/package=terra}
#' @seealso \code{\link{read_ecoreg}}() and \code{\link{make_ecoreg}}() to
#' prepare the ecoregion layer used here; \code{\link{cv_range}}() and
#' \code{\link{evaluate_range}}() to evaluate the resulting range map;
#' \code{\link{merge_range}}() to dissolve the returned polygons into a
#' single range polygon.
#' @example inst/examples/get_range_help.R
#' @importFrom methods is new
#' @importFrom terra vect crds intersect buffer makeValid aggregate rast disagg rasterize crop
#' @importFrom FNN knn.dist
#' @importFrom stats kmeans
#' @importFrom mclust Mclust mclustBIC emControl
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
                       format = c("SpatVector", "sf", "SpatRaster"),
                       res = 0.1,
                       verbose = TRUE) {

  ######################################################
  ### Stop messages
  ######################################################


  check_numeric(degrees_outlier, "degrees_outlier")
  check_numeric(clust_pts_outlier, "clust_pts_outlier")
  check_numeric(buff_width_point, "buff_width_point")
  check_numeric(buff_incrmt_pts_line, "buff_incrmt_pts_line")
  check_numeric(buff_width_polygon, "buff_width_polygon")
  format <- match.arg(format)
  check_numeric(res, "res")
  check_logical(verbose, "verbose")

  # Convex hulls are written to a PRIVATE sub-directory of 'dir_temp'. The
  # caller's directory is never deleted: 'dir_temp' defaults to tempdir(),
  # and removing that would destroy the session's temporary files. on.exit()
  # guarantees cleanup even if the function exits with an error.
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)
  }
  run_dir <- file.path(dir_temp, basename(tempfile("gbif_range_")))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

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
    message(
      paste0(
        "## Start of computation for species: ", sp.name, " ###\n"
      ),
      appendLF = FALSE
    )
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
    message(
      paste0(sum(cond)," outlier's from " ,
        nrow(occ_coord), " | proportion from total points: ",
        round((sum(cond) / nrow(occ_coord)) * 100,0), "%"
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

  # Restrict the layer to the occupied ecoregions and repair geometry once.
  # Simplification used to happen here (and before that, once per loop
  # iteration). It was dropped: on a global ecoregion layer it removed only
  # ~5% of vertices for ~6% of total runtime, and because Douglas-Peucker
  # works feature by feature it displaced each shared border independently of
  # its neighbour's, which is what produced slivers between adjacent
  # ecoregions. makeValid() keeps the validity repair without moving vertices.
  eco_all <- as.data.frame(ecoreg)[, ecoreg_name]
  keep_eco <- eco_all %in% uniq
  keep_eco[is.na(keep_eco)] <- FALSE

  ecoreg_s <- terra::makeValid(ecoreg[keep_eco, ])

  # Attribute column used to split ecoregions (hoisted: as.data.frame() on the
  # full layer inside the loop rebuilt the whole table on every iteration)
  eco_vals <- as.data.frame(ecoreg_s)[, ecoreg_name]
  
  # Loop over ecoregions
  SP.dist <- list()
  for (g in seq_along(uniq)) {
    
    # Print
    if (verbose){
      message(
        paste0(
          "ecoregion ", g, " of ", length(uniq), ": ", uniq[g], "\n"
        ),
        appendLF = FALSE
      )
    }

    # Handle NAs
    q1 <- eco_vals == uniq[g]
    q1[is.na(q1)] <- FALSE

    # Ecoregion polygon (already simplified and made valid before the loop)
    tmp <- ecoreg_s[q1, ]

    # Zero-width buffer once per ecoregion, not once per cluster
    b2 <- terra::buffer(tmp, width = 0)

    a <- ovo.coord.mod[ovo.coord.mod[[ecoreg_name]][[1]] == uniq[g]]

    # Coordinates are reused several times below; extract once
    xy <- terra::crds(a)

    # Jitter degenerate configurations so clustering can run. The previous
    # test (all longitudes equal to their own latitude) never fired, and its
    # body assigned to terra::crds(a)[, 2], for which no replacement method
    # exists. The real degenerate case is too few distinct locations, and the
    # offset must vary by row to actually separate coincident points. It is
    # deterministic so that results stay reproducible.
    if (nrow(unique(xy)) < 3) {
      xy[, 2] <- xy[, 2] + seq_len(nrow(xy)) * 0.00001
    }
    
    if (length(a) < 3) {
      k <- 1
      cluster.k <- stats::kmeans(
        x = xy,
        centers = k
      )
      cluster.k$clusters <- cluster.k$cluster 

    } else {

      # Determine number of clusters. mclust's default iteration limit is
      # effectively unlimited (.Machine$integer.max for both the EM loop and
      # the inner M-step loop used by VEI/VEE/VVE/VEV), so a non-converging
      # fit can run indefinitely. Capping both guarantees termination.
      m.clust <- mclust::Mclust(
        data = xy + 1000,
        verbose = FALSE,
        control = mclust::emControl(itmax = c(1000, 100))
      )
      
      # k = number of clusters (fall back to 1 if no model could be fitted)
      k <- if (is.null(m.clust)) 1L else m.clust$G
      
      # Reduce k if necessary so that KMeans_rcpp() will run
      while (k > length(a)-2) {k <- k-1} 
      if (k == 0) {k <- 1}
      
      cluster.k <- ClusterR::KMeans_rcpp(
        data = xy,
        clusters = k,
        num_init = 20,
        initializer = 'random'
      )
      
      while (length(unique(cluster.k$clusters)) < k) {
        k <- k-1
        cluster.k <- ClusterR::KMeans_rcpp(
          data = xy,
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
        temp_dir = run_dir,
        g = g,
        verbose = verbose
      )
      
      # Intersect polygon with ecoregion (zero buffer to avoid error).
      # 'b2' is constant within this ecoregion and is computed above.
      b1 <- terra::buffer(my.shpe, width = 0)
      polygons.list[[i]] <- terra::intersect(b1, b2)
    }  
    
    SP.dist[[g]] <- do.call("rbind", polygons.list)
  } 
  
  # Drop ecoregions that produced nothing. is.na() on a list of SpatVectors
  # is not a reliable emptiness test.
  keep <- vapply(
    SP.dist,
    function(z) !is.null(z) && nrow(z) > 0,
    logical(1)
  )
  lala <- SP.dist[keep]
  

  ### =========================================================================
  ### Check and return output
  ### =========================================================================
  

  if (length(lala) == 0) {
    stop('No occurrences within ecoregions. Empty raster produced...')
  }
  shp.species <- do.call("rbind", lala)

  # Keep only the ecoregion label and the species name, so that the returned
  # schema does not depend on which columns the user's 'ecoreg' layer carried.
  # Guarded, because terra::intersect() renames fields on name collisions.
  if (ecoreg_name %in% names(shp.species)) {
    shp.species <- shp.species[, ecoreg_name]
  }
  # Remove the overlaps between buffered hulls of neighbouring clusters by
  # dissolving WITHIN each ecoregion. All pieces of a given ecoregion were
  # clipped against the same geometry, so this union has no mismatched borders.
  # Across ecoregions nothing is merged, which is where slivers would arise.
  # terra::erase() was used here previously, but the SpatVector it returns is
  # internally corrupted (terra issue #1710): it prints and validates, yet
  # fails on column assignment, rasterize() and writeVector().
  shp.species <- terra::makeValid(shp.species)
  shp.species <- terra::aggregate(
    x = shp.species,
    by = ecoreg_name,
    count = FALSE
  )
  shp.species$species <- occ_coord$input_search[1]

  # Convert to requested format
  if (format == "SpatRaster") {
    res.use <- 1 / res
    ras.res <- terra::rast(terra::disagg(terra::rast(), res.use))
    ras <- terra::rasterize(shp.species, ras.res)
    shp.species <- terra::crop(ras, shp.species)
    names(shp.species) <- occ_coord$input_search[1]

  } else if (format == "sf") {
    shp.species <- sf::st_as_sf(shp.species)
  }
  # SpatVector: returned as-is, attributes already set
  
  # Final print
  if (verbose){
    message(
      paste0(
        "## End of computation for species: ", sp.name, " ###\n"
      ),
      appendLF = FALSE
    )
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
                           format = format,
                           res = res)
  result$rangeOutput <- shp.species

  return(result)
}
