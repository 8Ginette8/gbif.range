### =========================================================================
### make_tiles
### =========================================================================
#' Generate Spatial tiles from a raster extent or geometry
#'
#' Divides a given spatial extent or geometry into a specified number of tiles.
#' Each tile is returned as a \code{POLYGON()} geometry. The original extent
#' can be: (1) converted into a single \code{POLYGON()}, or (2) subdivided
#' into approximately \code{ntiles} regular fragments, each returned as a
#' \code{POLYGON()} and, optionally, as a \code{SpatExtent}. Based on a
#' specific extent, one or several tiles are generated. Tiles can be smaller
#' raster extents or geometry arguments \code{POLYGON()}. The original extent
#' is therefore either converted into a \code{POLYGON()} argument, or divided
#' into \code{ntiles} of regular fragments which are converted into
#' \code{POLYGON()} arguments and smaller SpatExtent.
#'
#' @param geo Object of class \code{Extent}, \code{SpatExtent},
#' \code{SpatialPolygon}, \code{SpatialPolygonDataframe} or \code{SpaVector}
#' (WGS84 or planar) to define the study's area extent.
#' Default is \code{NULL} i.e. the whole globe.
#' @param ntiles Numeric. In how many tiles/fragments should geo be
#' divided approximately?
#' @param sext Logical. Should a list of SpatExtent also be returned
#' for each generated \code{POLYGON()}?
#' @return A list of geometry arguments POLYGON() of length \code{ntiles},
#' and of \code{SpatExtent} if \code{sext = TRUE}).
#' @references 
#' Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P.,
#' Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate,
#' soil, and land cover on plant species distribution in the European Alps.
#' Ecological monographs, 91(2), e01433. 10.1002/ecm.1433
#' @example inst/examples/make_tiles_help.R
#' @importFrom terra ext
#' @export
make_tiles <- function(geo, ntiles, sext = TRUE){
  
    ###########################################
    ### Check input data
    ###########################################
t

    # For study area
	if (!is.null(geo)) {
		if (!(class(geo) %in% "SpatExtent")) {geo <- terra::ext(geo)}
	} else {
		geo <- terra::ext(-180, 180, -90, 90)
	}

    # General
    check_numeric(ntiles, "ntiles")
    check_logical(sext, "sext")


    ###########################################
    ### Check input data
    ###########################################


	# Divide original extent into smaller ones otherwise
	n.tiles <- seq(1:sqrt(ntiles))
	dtiles <- data.frame(n.tiles)
	xFactor <- (geo$xmax - geo$xmin) / length(n.tiles)
	yFactor <- (geo$ymax - geo$ymin) / length(n.tiles)
	dtiles$xCH <- -dtiles$n.tiles * xFactor + geo$xmax
	dtiles$yCH <- -dtiles$n.tiles * yFactor + geo$ymax
	allxy <- rbind(c(geo$xmax, geo$ymax), dtiles[, -1])
	tile.index <- nrow(allxy) - 1

	# Create all smaller extents
	all.tiles <-
	lapply(1:tile.index, function(x){

		# Set target ymax and ymin
		t.ymax <- allxy[x, 2]
		t.ymin <- allxy[x + 1, 2]

		# Return NULL if t.ymax = t.ymin
		if (t.ymax == t.ymin) {return(list(NULL, NULL))}

		# Generate tiles for each line (xmax -> xmin)
		line.tile <-
		lapply(1:tile.index, function(y){

			# Set target xmax and xmin
			t.xmax <- allxy[y, 1]
			t.xmin <- allxy[y + 1, 1]

			# Return NULL if t.xmax = t.xmin
			if (t.xmax == t.xmin) {return(list(NULL, NULL))}

			# Generate the tile
			one.tile <- paste0("POLYGON((", t.xmin, " ", t.ymin, ", ",
				 						t.xmax, " ", t.ymin, ", ",
				 						t.xmax, " ", t.ymax, ", ",
			 							t.xmin, " ", t.ymax, ", ",
			 							t.xmin, " ", t.ymin, "))")
			meta.tile <- c(t.xmin, t.xmax, t.ymin, t.ymax)
			return(list(one.tile, meta.tile))
		})
		part.tile <- lapply(line.tile, function(y) y[[1]])
		part.meta <- lapply(line.tile, function(y) y[[2]])
		return(list(part.tile, part.meta))
	})

	# Unlist tile geo
	part.tile2 <- lapply(all.tiles, function(x) x[[1]])
	geo.tiles <- unlist(part.tile2, recursive = FALSE)
	geo.tiles[vapply(geo.tiles, is.null, logical(1))] <- NULL

	# Return if xmin == xmax
	if (is.null(unlist(geo.tiles))) {return(NULL)}

	# Unlist tile meta
	geo.meta <- unlist(lapply(all.tiles, function(x) x[[2]]), recursive = FALSE)
	geo.meta[vapply(geo.meta, is.null, logical(1))] <- NULL
	geo.meta <- lapply(geo.meta, function(x) terra::ext(x))

	# Return
	if (sext) {
		return(list(geo.tiles, geo.meta))
	} else {
		return(geo.tiles)
	}
}
