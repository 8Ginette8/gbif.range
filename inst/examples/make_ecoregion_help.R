# Open data
rst.path <- paste0(system.file(package = "gbif.range"), "/extdata/rst_enl.tif")
rst <- terra::rast(rst.path)
shp.path <- paste0(system.file(package = "gbif.range"), "/extdata/shp_lonlat.shp")
shp.lonlat <- terra::vect(shp.path)
rst <- terra::crop(rst, shp.lonlat)

# Apply the function by infering 50 classes of environments
my.eco <- make_ecoregion(env = rst,
						 nclass = 50)

# Test plot
terra::plot(my.eco)

# Downloading in the European Alps the observations of one plant species
obs.arcto <- get_gbif(sp_name = "Arctostaphylos alpinus",
					  geo = shp.lonlat)

# Create the range map based on our custom ecoregion
range.arcto <- get_range(occ_coord = obs.arcto,
						 bioreg = my.eco, 
						 bioreg_name = "EcoRegion",
						 res = 20)

# Plot
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(range.arcto$rangeOutput,add = TRUE,
	col = "darkgreen", axes = FALSE, legend = FALSE)
graphics::points(obs.arcto[, c("decimalLongitude","decimalLatitude")],
	pch = 20, col = "#99340470", cex= 1)