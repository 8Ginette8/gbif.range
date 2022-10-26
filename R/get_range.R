##################################################
## get_range
##################################################
#' Create a species range map based on a get_gbif() output
#' 
#' Estimates species ranges based on occurrence data (GBIF or not) and bioregions.
#' It first deletes outliers from the observation dataset and then creates a polygon
#' (convex hull) with a user specified buffer around all the observations of one bioregion.
#' If there there is only one observation in a bioregion, a buffer around this point
#' will be created. If all points in a bioregion are on a line, the function will also
#' create a buffer around these points, however, the buffer size increases with the number
#' of points in the line.
#' 
#' @param species_name character string of the species name. E.g. "Anemone nemorosa".
#' @param occ_coord a SpatialPoints object.
#' @param Bioreg shapefile containg different bioregions (convex hulls will be classified
#' on a bioreg basis).
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
#' @details ...
#' @return a shapefile.
#' @references 
#' Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J. S., ... &
#' Pellissier, L. (2022). An integrated high‐resolution mapping shows congruent biodiversity
#' patterns of Fagales and Pinales. New Phytologist, 235(2), 759-772 10.1111/nph.18158
#' 
#' Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Terra - CRAN
#' @seealso ...
#' @examples
#' @export
get_range <- function (species_name = NULL, 
                       occ_coord = NULL, 
                       Bioreg = NULL, 
                       Bioreg_name = NULL, 
                       degrees_outlier=10,
                       clustered_points_outlier=3,
                       buffer_width_point=0.5, 
                       buffer_increment_point_line=0.5, 
                       buffer_width_polygon=0.1, 
                       dir_temp=paste0("temp",sample(1:99999999,1))){
  
  ### =========================================================================
  ### remove duplicates
  ### =========================================================================

  occ_coord@coords <- round(occ_coord@coords, 4)
  occ_coord <- remove.duplicates(occ_coord)
  
  ### =========================================================================
  ### Check if sufficient data
  ### =========================================================================
  
  # Check if there sufficient species & if not, make an entry in the log-file and
  # end the function
  if (length(occ_coord)<=clustered_points_outlier+1){
    stop("Too few occurences!")
  } 
    
  cat("## Start of computation for species: ",species_name," ###", "\n") 
  
  ### =========================================================================
  ### Identify outliers
  ### =========================================================================
  
  #create distance matrix...
  mat_dist <- as.matrix(knn.dist(occ_coord@coords, k=clustered_points_outlier))
  
  #mark outliers
  cond <- apply(mat_dist, 1, function(x) x[clustered_points_outlier])>degrees_outlier
  rm(mat_dist) 
  
  cat(paste0(sum(cond), " outlier's from " ,nrow(occ_coord), " | proportion from total points: ", round((sum(cond)/length(occ_coord))*100,0), "%\n"))
  
  occ_coord_mod <- occ_coord[!cond,]

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
  
  cat("## End of computation for species: ",species_name," ###", "\n") 

  return(shp_species)
}