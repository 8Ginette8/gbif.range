# Open data
rst.path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/rst_enl.tif"
)
rst <- terra::rast(rst.path)
shp.path <- paste0(
	system.file(package = "gbif.range"),
	"/extdata/shp_lonlat.shp"
)
shp.lonlat <- terra::vect(shp.path)
rst <- terra::crop(rst, shp.lonlat)

# Apply the function by inferring 200 environmental classes
my.eco <- make_ecoreg(env = rst,
	nclass = 200,
	format = "sf"
)

\donttest{
# Downloading in the European Alps the observations of one plant species
obs.paed <- get_gbif(
  sp_name = "Paederota bonarota",
  geo = shp.lonlat,
  grain = 1
)

# Create the range map based on:
# - custom ecoregion at 5 x 5 km resolution
# - smaller buffer because of regional extent
range.arcto <- get_range(
	occ_coord = obs.paed,
	ecoreg = my.eco,
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
alps.shp <- terra::aggregate(terra::crop(countries,terra::ext(rst)))
r.arcto <- terra::mask(range.arcto$rangeOutput,alps.shp)
terra::plot(terra::crop(countries,terra::ext(rst)), col = "#bcbddc")
terra::plot(
	r.arcto,
	add = TRUE,
	col = "darkgreen",
	axes = FALSE,
	legend = FALSE
)
graphics::points(
	obs.paed[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#99340470",
	cex = 1
)

}
