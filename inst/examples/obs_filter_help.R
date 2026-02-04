# Load data
shp.path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/shp_lonlat.shp"
)
shp.lonlat <- terra::vect(shp.path)
rst.path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/rst_enl.tif"
)
rst <- terra::rast(rst.path)

# Downloading in the European Alps the observations of two plant species
obs.arcto <- get_gbif(
	sp_name = "Arctostaphylos alpinus",
	geo = shp.lonlat
)
obs.saxi <- get_gbif(
	sp_name = "Saxifraga cernua",
	geo = shp.lonlat
)

# Test plot
terra::plot(shp.lonlat)
graphics::points(
	obs.arcto[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#238b4550",
	cex = 1
)
graphics::points(
	obs.saxi[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#99000d50",
	cex = 1
)

# rbind both datasets
both.sp <- rbind(obs.arcto,obs.saxi)

# Run function
obs.filt <- obs_filter(gbifs = both.sp, grid = rst, threshold = 4)

# Check new points
terra::plot(shp.lonlat)
graphics::points(
	obs.filt[obs.filt$Species%in%"Arctostaphylos alpinus",c("x","y")],
	pch = 20,
	col = "#238b4550",
	cex = 1
)
graphics::points(
	obs.filt[obs.filt$Species%in%"Saxifraga cernua",c("x","y")],
	pch = 20,
	col = "#99000d50",
	cex = 1
)