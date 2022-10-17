### =========================================================================
### make_tiles
### =========================================================================
#' Create a specific number of tiles based on an extent
#'
#' Not to be called directly by the user
#'
#' @author Yohann Chauvier
#' @export
#' 
make_tiles = function(geo, Ntiles, meta = TRUE){

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

	# Return
	if (meta) {
		return(list(geo.tiles,geo.meta))
	} else {
		return(geo.tiles)
	}
}

