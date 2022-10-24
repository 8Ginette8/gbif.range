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
#' @param geo Object of class 'Extent', 'SpatExtent', 'SpatialPolygon', 'SpatialPolygonDataframe',
#' or 'SpaVector' (WGS84 or planar) to define the study's area extent. Default is NULL i.e. the
#' whole globe.
#' @param Ntiles Numeric. In how many tiles/fragments should geo be divided approximately?
#' @param sext Logical. Should a list of SpatExtent also be returned for each generated POLYGON()?
#' @return A list of geometry arguments POLYGON() of length Ntiles (and of SpatExtent
#' if sext=TRUE)
#' @references 
#' Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann,
#' N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the
#' European Alps. Ecological monographs, 91(2), e01433. 10.1002/ecm.1433
#' @examples
#' 
#' # Load the European Alps Extent
#' data(geo_dat)
#' 
#' # Apply the function to divide the extent in ~20 fragments
#' mt = make_tiles(geo=shp.lonlat,Ntiles=20,sext=TRUE); mt
#' @export
#' 
make_tiles = function(geo, Ntiles, sext = TRUE){

	# For study area
	if (!is.null(geo)) {
		if (!(class(geo)%in%"SpatExtent")) {geo = ext(geo)}
	} else {
		geo = ext(-180,180,-90,90)
	}

	# Divide original extent into smaller ones otherwise
	ntiles = seq(1:sqrt(Ntiles))
	dtiles = data.frame(ntiles)
	xFactor = (geo$xmax - geo$xmin)/length(ntiles)
	yFactor = (geo$ymax - geo$ymin)/length(ntiles)
	dtiles$xCH = -dtiles$ntiles*xFactor + geo$xmax
	dtiles$yCH = -dtiles$ntiles*yFactor + geo$ymax
	allxy = rbind(c(geo$xmax,geo$ymax),dtiles[,-1])
	tile.index = nrow(allxy)-1

	# Create all smaller extents
	all.tiles =
	lapply(1:tile.index, function(x){

		# Set target ymax and ymin
		t.ymax = allxy[x,2]
		t.ymin = allxy[x+1,2]

		# Generate tiles for each line (xmax -> xmin)
		line.tile =
		lapply(1:tile.index,function(y){

			# Set target xmax and xmin
			t.xmax = allxy[y,1]
			t.xmin = allxy[y+1,1]

			# Generate the tile
			one.tile = paste0("POLYGON((",t.xmin," ",t.ymin,", ",
				 						t.xmax," ",t.ymin,", ",
				 						t.xmax," ",t.ymax,", ",
			 							t.xmin," ",t.ymax,", ",
			 							t.xmin," ",t.ymin,"))")
			meta.tile = c(t.xmin,t.xmax,t.ymin,t.ymax)
			return(list(one.tile,meta.tile))
		})
		part.tile = lapply(line.tile,function(y) y[[1]])
		part.meta = lapply(line.tile,function(y) y[[2]])
		return(list(part.tile,part.meta))
	})

	# Unlist tile geo
	part.tile2 = lapply(all.tiles,function(x) x[[1]])
	geo.tiles = unlist(part.tile2,recursive=FALSE)

	# Unlist tile meta
	geo.meta = unlist(lapply(all.tiles,function(x) x[[2]]),recursive=FALSE)
	geo.meta = lapply(geo.meta,function(x) ext(x))

	# Return
	if (sext) {
		return(list(geo.tiles,geo.meta))
	} else {
		return(geo.tiles)
	}
}

