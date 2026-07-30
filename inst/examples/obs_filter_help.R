\donttest{
# Load data
shp_path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/shp_lonlat.shp"
)
shp_lonlat <- terra::vect(shp_path)
rst_path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/rst_enl.tif"
)
rst <- terra::rast(rst_path)

# Download observations for two plant species in the European Alps
obs_paed <- get_gbif(
	sp_name = "Paederota bonarota",
	geo = shp_lonlat
)
obs_saxi <- get_gbif(
	sp_name = "Saxifraga cernua",
	geo = shp_lonlat
)

# Guard against an unavailable remote source
if (gbif_have(obs_paed, obs_saxi)) {

# Test plot
terra::plot(shp_lonlat)
graphics::points(
	obs_paed[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#238b4550",
	cex = 1
)
graphics::points(
	obs_saxi[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#99000d50",
	cex = 1
)

# Combine both datasets
both_sp <- rbind(obs_paed, obs_saxi)

# Run function
obs_filt <- obs_filter(gbifs = both_sp, grid = rst, threshold = 4)

# Check new points
terra::plot(shp_lonlat)
graphics::points(
	obs_filt[obs_filt$Species%in%"Paederota bonarota", c("x","y")],
	pch = 20,
	col = "#238b4550",
	cex = 1
)
graphics::points(
	obs_filt[obs_filt$Species%in%"Saxifraga cernua", c("x","y")],
	pch = 20,
	col = "#99000d50",
	cex = 1
)

}
}
