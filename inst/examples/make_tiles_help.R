# Load the European Alps Extent
shp_path <- paste0(
	system.file(package = "gbif.range"), 
	"/extdata/shp_lonlat.shp"
)
shp_lonlat <- terra::vect(shp_path)
 
# Apply the function to divide the extent in ~20 fragments
mt <- make_tiles(geo = shp_lonlat, ntiles = 20, sext = TRUE)
