# Load the European Alps Extent
shp.path <- paste0(system.file(package = "gbif.range"), "/extdata/shp_lonlat.shp")
shp.lonlat <- terra::vect(shp.path)
 
# Apply the function to divide the extent in ~20 fragments
mt <- make_tiles(geo = shp.lonlat,
				 Ntiles = 20,
				 sext = TRUE)