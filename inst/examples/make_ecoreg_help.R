# Open data
rst_path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/rst_enl.tif"
)
rst <- terra::rast(rst_path)
shp_path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/shp_lonlat.shp"
)
shp_lonlat <- terra::vect(shp_path)
rst <- terra::crop(rst, shp_lonlat)

# Apply the function by inferring 200 environmental classes
my_eco <- make_ecoreg(env = rst,
	nclass = 200,
	format = "sf"
)

\donttest{
# Downloading in the European Alps the observations of one plant species
obs_paed <- get_gbif(
  sp_name = "Paederota bonarota",
  geo = shp_lonlat,
  grain = 1
)

# Guard against an unavailable remote source
if (gbif_have(obs_paed)) {

# Create the range map based on:
# - custom ecoregion at 5 x 5 km resolution
# - smaller buffer because of regional extent
range_paed <- get_range(
	occ_coord = obs_paed,
	ecoreg = my_eco,
	ecoreg_name = "EcoRegion",
	res = 0.05,
  degrees_outlier = 0.5,
  buff_width_point = 0.5,
  buff_incrmt_pts_line = 0.5,
  buff_width_polygon = 0.5
)

# Plot
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(
	merge_range(range_paed),
	add = TRUE,
	col = "darkgreen",
	axes = FALSE,
	legend = FALSE
)
graphics::points(
	obs_paed[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#99340470",
	cex = 1
)

}
}
