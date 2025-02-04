# Open
rst_path <- paste0(system.file(package = "gbif.range"), "/extdata/rst_enl.tif")
rst <- terra::rast(rst_path)
shp_path <- paste0(system.file(package = "gbif.range"), "/extdata/shp_lonlat.shp")
shp_lonlat <- terra::vect(shp_path)
rst <- terra::crop(rst, shp_lonlat)

# Apply the function by infering 50 classes of environments
my_eco <- make_ecoregion(rst, nclass = 50)
terra::plot(my_eco)

# Downloading in the European Alps the observations of one plant species
obs_arcto <- get_gbif("Arctostaphylos alpinus", geo = shp_lonlat)

# Create the range map based on our custom ecoregion
range_arcto <- get_range(obs_arcto, my_eco, "EcoRegion", res = 20)

# Plot
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(range_arcto$range_output,add = TRUE,
	col = "darkgreen", axes = FALSE, legend = FALSE)
graphics::points(obs_arcto[, c("decimalLongitude","decimalLatitude")],
	pch = 20, col = "#99340470", cex= 1)