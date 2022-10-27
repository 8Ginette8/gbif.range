##################################################
## get_range
##################################################
#' Create a species range map based on a *get_gbif()* output
#' 
#' Estimates species ranges based on occurrence data (GBIF or not) and bioregions.
#' It first deletes outliers from the observation dataset and then creates a polygon
#' (convex hull) with a user specified buffer around all the observations of one bioregion.
#' If there is only one observation in a bioregion, a buffer around this point
#' will be created. If all points in a bioregion are on a line, the function will also
#' create a buffer around these points, however, the buffer size increases with the number
#' of points in the line.
#' 
#' @param sp_name Character of the species name. E.g. "Anemone nemorosa".
#' @param occ_coord *get_gbif()* output or a SpatVector/SpatialPoints object.
#' @param Bioreg A SpatVector or SpatialPolygonsDataFrame containg different bioregions
#' (convex hulls will be classified on a bioreg basis). Although whatever shapefile may be
#' set as input, note that three ecoregion shapefiles are already included in the library:
#' *eco.earh* (for terrestrial species; Nature conservancy version adapted from Olson & al. 2001),
#' *eco.marine* (for coastal and reef species; Spalding & al. 2007) and *eco.fresh*
#' (for freshwater species; Abell & al. 2008). For deep ocean/sea species, *eco.earth* may be
#' used, but the polygon estimates will only be geographic.
#' @param Bioreg_name how is the slot containing the bioregion names called?
#' @param degrees_outlier distance threshold (degrees) for outlier classification. If the
#' nearest minimal distance to the next point is larger than this threshold, it will be
#' considered as an outlier.
#' @param clustered_points_outlier maximum number of points which are closer to each other
#' than the degrees_outlier, but should still be considered as outliers.
#' @param buffer_width_point buffer (in degrees) which will be applied around single observations.
#' @param buffer_increment_point_line how much should the buffer be increased for each point on a line.
#' @param buffer_width_polygon buffer (in degrees) which will be applied around distribution polygons
#' (for each bioregion).
#' @param dir_temp where should the temporary text file for the convex hull be saved?
#' (text file will be deleted again).
#' @param raster Logical. Should the output be a unified raster? Default is TRUE
#' @param res Numeric. If raster = TRUE, which resolution? Final resolution in ° = 1°/res
#' e.g.,  = 0.1° if res = 10. Default is 10.
#' @details ...
#' @return A shapefile or a SpatRaster
#' @references
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
#' Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A. Ferdaña, Max
#' Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana, Sara A. Lourie, Kirsten D.
#' Martin, Edmund McManus, Jennifer Molnar, Cheri A. Recchia, James Robertson, Marine
#' Ecoregions of the World: A Bioregionalization of Coastal and Shelf Areas, BioScience,
#' Volume 57, Issue 7, July 2007, Pages 573–583. doi: 10.1641/B570707
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
#' @seealso ...
#' @examples
#' # Load available ecoregions
#' data(ecoregions)
#' 
#' # First download the worldwide observations of Panthera tigris and convert to SpatialPoints
#' obs.pt = get_gbif("Panthera tigris",basis=c("OBSERVATION","HUMAN_OBSERVATION"))
#' sp.shp = SpatialPoints(obs.pt[c("decimalLongitude","decimalLatitude")],proj4string=crs(eco.earth))
#' 
#' # Plot
#' plot(eco.earth)
#' plot(sp.shp,pch=20,col="#238b4550",cex=4,add=TRUE)
#' 
#' # Generate the distributional range map of Panthera tigris for the finest terrestrial ecoregions
#' range.tiger = get_range("Panthera tigris",sp.shp,eco.earth,"ECO_NAME")
#' 
#' # Plotting
#' plot(eco.earth)
#' plot(range.tiger,col="#238b45",add=TRUE)
#' 
#' @export
get_range <- function (sp_name = NULL, 
                       occ_coord = NULL, 
                       Bioreg = NULL, 
                       Bioreg_name = NULL, 
                       degrees_outlier = 3,
                       clustered_points_outlier = 2,
                       buffer_width_point = 4, 
                       buffer_increment_point_line = 0.5, 
                       buffer_width_polygon = 4,
                       dir_temp = paste0("temp",sample(1:99999999,1)),
                       raster = TRUE,
                       res = 10){

  ### =========================================================================
  ### Convert all possible object in SpatVector
  ### =========================================================================

  if (class(occ_coord)%in%"SpatialPoints") {
    occ_coord <- vect(occ_coord)
  } else if (any(names(occ_coord)%in%"decimalLongitude")) {
    occ_coord <- vect(occ_coord[,c("decimalLongitude","decimalLatitude")],
      geom=c("decimalLongitude","decimalLatitude"),crs="+init=epsg:4326")
  }

  ### =========================================================================
  ### remove duplicates
  ### =========================================================================

  occ_coordS <- crds(occ_coord)
  occ_coordS <- round(occ_coordS, 4)
  occ_coordS <- occ_coordS[!duplicated(occ_coords),]
  
  ### =========================================================================
  ### Check if sufficient data
  ### =========================================================================
  
  # Check if there sufficient species & if not, make an entry in the log-file and
  # end the function
  if (nrow(occ_coordS)<=clustered_points_outlier+1){
    stop("Too few occurences!")
  } 
    
  cat("## Start of computation for species: ",sp_name," ###", "\n") 
  
  ### =========================================================================
  ### Identify outliers
  ### =========================================================================
  
  #create distance matrix...
  mat_dist <- as.matrix(knn.dist(occ_coordS, k=clustered_points_outlier))
  
  #mark outliers
  cond <- apply(mat_dist, 1, function(x) x[clustered_points_outlier])>degrees_outlier
  rm(mat_dist) 
  
  cat(paste0(sum(cond), " outlier's from " ,nrow(occ_coordS), " | proportion from total points: ", round((sum(cond)/length(occ_coord))*100,0), "%\n"))
  
  occ_coord_mod <- occ_coordS[!cond,]

  # Stop if too many outliers in data set
  if(length(occ_coord_mod)==0){
    stop('Too few occurrences within outlier threshold!')
  } 

  ### =========================================================================
  ### Define those polygons
  ### =========================================================================
  
  ovo=over(occ_coord_mod,Bioreg)
  uniq <- levels(factor(ovo[,Bioreg_name]))
  
  # loop over bioregions
  SP_dist <- list()
  for(g in 1:length(uniq)) {
    
    cat('Bioregion', g, ' of ',length(uniq),": ",uniq[g], '\n')
    
    tmp <- as(gSimplify(Bioreg[Bioreg@data[[Bioreg_name]] == uniq[g],],tol=0.001,topologyPreserve=TRUE),"SpatialPolygons")
    a=occ_coord_mod[which(ovo[,Bioreg_name]==uniq[g])]
    
    if(length(a)<3){
      k <- 1
      cluster_k <- kmeans(a@coords,k)
      cluster_k$clusters <- cluster_k$cluster 
    } else {
      if(all(a@coords[,1]==a@coords[,2])){
        a@coords[,2] <- a@coords[,2]+0.00001
      }
      
      m_clust <- Mclust(a@coords+1000, verbose=F) #to determine number of clusters
      k <- m_clust$G #k=number of clusters
      
      while(k>length(a)-2){k <- k-1} #reduce k if necessary so that KMeans_rcpp() will run
      if(k==0){k <- 1}
      
      cluster_k <- KMeans_rcpp(a@coords, k, num_init = 20, initializer = 'random')
      
      while(length(unique(cluster_k$clusters))<k){
        k <- k-1
        cluster_k <- KMeans_rcpp(a@coords, k, num_init = 20, initializer = 'random')
      }
      
    }
   
    polygons_list <- list() 
    for(i in 1:k){
      a_temp <- a[cluster_k$clusters==i] #kmeans (with number of clusters from mcluster)
      
      
      my_shpe=conv_function(a_temp,
                            bwp=buffer_width_point,
                            bipl=buffer_increment_point_line,
                            bwpo=buffer_width_polygon,
                            temp_dir=dir_temp,
                            g=g)
      
      polygons_list[[i]] <- suppressWarnings(gIntersection(gBuffer(SpatialPolygons(Srl=list(my_shpe)), width=0),tmp)) #zero buffer to avoid error
      polygons_list[[i]]$ID <- i
    }  
    
    SP_dist[[g]] <- Reduce(rbind, polygons_list)

    if(class(SP_dist[[g]])=='SpatialCollections'){
      SP_dist[[g]] <- SP_dist[[g]]@polyobj #only keep SpatialPolygon
    }
    
  } 
  
  L <- SP_dist[!is.na(SP_dist)]
  
  ### =========================================================================
  ### Check and return output
  ### =========================================================================
  
  if(!dir.exists(dir_temp)){
    unlink(dir_temp, recursive=T)
  }
  
  if(length(L)==0){
    stop('No occurrences within Bioregions. Empty raster produced.')
  } 
    
  shp_species <- Reduce(rbind, L)
  # shp_species=gUnaryUnion(shp_species,checkValidity = 2)
  shp_species@proj4string=occ_coord@proj4string

  cat("## End of computation for species: ",sp_name," ###", "\n") 

  # Convert in raster files or not
  if (raster){
    ras.res <- rast(disaggregate(raster(),res))
    sp.range.u <- gUnaryUnion(shp_species)
    ras <- rasterize(vect(sp.range.u),ras.res)
    shp_species <- crop(ras,ext(sp.range.u))
  }
  
  cat("## End of computation for species: ",sp_name," ###", "\n") 

  return(shp_species)
}