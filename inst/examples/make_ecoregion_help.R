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
col_palette = colorRampPalette(c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
            "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
            "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"))
colcol = col_palette(length(my.eco))
my.eco$color = sample(paste0(colcol, ""), length(my.eco), replace = FALSE)
terra::plot(my.eco, col = my.eco$color)

# Downloading in the European Alps the observations of one plant species
obs.arcto <- get_gbif(sp_name = "Arctostaphylos alpinus",
		      geo = shp.lonlat)

# Create the range map based on our custom ecoregion at 5 x 5 km resolution
range.arcto <- get_range(occ_coord = obs.arcto,
			 bioreg = my.eco, 
			 bioreg_name = "EcoRegion",
			 res = 0.05)

# Plot
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(range.arcto$rangeOutput,add = TRUE,
	col = "darkgreen", axes = FALSE, legend = FALSE)
graphics::points(obs.arcto[, c("decimalLongitude","decimalLatitude")],
	pch = 20, col = "#99340470", cex= 1)
