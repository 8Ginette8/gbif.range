### =========================================================================
### make_tiles
### =========================================================================
#' Create Tiled GBIF Geometry Queries
#'
#' Divide a study extent into a set of smaller \code{POLYGON()} query strings
#' that can be used with the GBIF \code{geometry} parameter.
#'
#' @param geo Spatial extent or geometry used to define the study area.
#' Accepted classes are \code{Extent}, \code{SpatExtent},
#' \code{SpatialPolygon}, \code{SpatialPolygonDataFrame}, \code{SpatVector},
#' and \code{sf}. If \code{NULL}, the whole globe is used.
#' @param ntiles Numeric approximate number of tiles to create.
#' @param sext Logical. Should the corresponding \code{SpatExtent} objects also
#' be returned?
#' @details The input extent is converted into a regular grid of smaller query
#' polygons. This is mainly intended to support tiled GBIF downloads when large
#' extents would otherwise return too many records in a single query.
#' @return If \code{sext = TRUE}, a list containing WKT \code{POLYGON()}
#' strings and matching \code{SpatExtent} objects. Otherwise, a vector of WKT
#' \code{POLYGON()} strings.
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

    # For study area
	if (!is.null(geo)) {
		if (inherits(geo, "sf")) geo <- terra::vect(geo)
		if (!inherits(geo, "SpatExtent")) geo <- terra::ext(geo)
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
