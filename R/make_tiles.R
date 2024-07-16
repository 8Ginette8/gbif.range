### =========================================================================
### make_tiles
### =========================================================================
#' Create a specific number of tiles based on a raster extent
#'
#' Based on a specific extent, one or several tiles are generated. Tiles can be smaller
#' raster extents or geometry arguments POLYGON(). The original extent is therefore either
#' converted into a POLYGON() argument, or divided into Ntiles of regular fragments which are
#' converted into POLYGON() arguments and smaller SpatExtent.
#'
#' @param geo Object of class Extent, SpatExtent, SpatialPolygon, SpatialPolygonDataframe,
#' or SpaVector (WGS84 or planar) to define the study's area extent. Default is NULL i.e. the
#' whole globe.
#' @param Ntiles Numeric. In how many tiles/fragments should geo be divided approximately?
#' @param sext Logical. Should a list of SpatExtent also be returned for each generated POLYGON()?
#' @return A list of geometry arguments POLYGON() of length Ntiles (and of SpatExtent
#' if sext=TRUE).
#' @references 
#' Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann,
#' N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the
#' European Alps. Ecological monographs, 91(2), e01433. 10.1002/ecm.1433
#' @example inst/examples/make_tiles_help.R
#' @importFrom terra ext
#' @export
#' 
make_tiles = function(geo, Ntiles, sext = TRUE){

	# For study area
	if (!is.null(geo)) {
		if (!(class(geo)%in%"SpatExtent")) {geo = terra::ext(geo)}
	} else {
		geo = terra::ext(-180,180,-90,90)
	}

	# Divide original extent into smaller ones otherwise
	ntiles = seq(1:sqrt(Ntiles))
	dtiles = data.frame(ntiles)
	xFactor = (geo$xmax - geo$xmin)/length(ntiles)
	yFactor = (geo$ymax - geo$ymin)/length(ntiles)
	dtiles$xCH = -dtiles$ntiles*xFactor + geo$xmax
	dtiles$yCH = -dtiles$ntiles*yFactor + geo$ymax
	allxy = rbind(c(geo$xmax,geo$ymax),dtiles[,-1])
	tile_index = nrow(allxy)-1

	# Create all smaller extents
	all_tiles =
	lapply(1:tile_index, function(x){

		# Set target ymax and ymin
		t_ymax = allxy[x,2]
		t_ymin = allxy[x+1,2]

		# Return NULL if t_ymax = t.ymin
		if (t_ymax == t_ymin) {return(list(NULL,NULL))}

		# Generate tiles for each line (xmax -> xmin)
		line_tile =
		lapply(1:tile_index,function(y){

			# Set target xmax and xmin
			t_xmax = allxy[y,1]
			t_xmin = allxy[y+1,1]

			# Return NULL if t.xmax = t.xmin
			if (t_xmax == t_xmin) {return(list(NULL,NULL))}

			# Generate the tile
			one_tile = paste0("POLYGON((",t_xmin," ",t_ymin,", ",
				 						t_xmax," ",t_ymin,", ",
				 						t_xmax," ",t_ymax,", ",
			 							t_xmin," ",t_ymax,", ",
			 							t_xmin," ",t_ymin,"))")
			meta_tile = c(t_xmin,t_xmax,t_ymin,t_ymax)
			return(list(one_tile,meta_tile))
		})
		part_tile = lapply(line_tile,function(y) y[[1]])
		part_meta = lapply(line_tile,function(y) y[[2]])
		return(list(part_tile,part_meta))
	})

	# Unlist tile geo
	part_tile2 = lapply(all_tiles,function(x) x[[1]])
	geo_tiles = unlist(part_tile2,recursive = FALSE)
	geo_tiles[sapply(geo_tiles, is.null)] = NULL

	# Return if xmin == xmax
	if (is.null(unlist(geo_tiles))) {return(NULL)}

	# Unlist tile meta
	geo_meta = unlist(lapply(all_tiles,function(x) x[[2]]),recursive = FALSE)
	geo_meta[sapply(geo_meta, is.null)] = NULL
	geo_meta = lapply(geo_meta,function(x) terra::ext(x))

	# Return
	if (sext) {
		return(list(geo_tiles,geo_meta))
	} else {
		return(geo_tiles)
	}
}