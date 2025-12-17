\dontrun{
###########################################
### Example plot
###########################################

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

# Apply the function by infering 50 classes of environments
my.eco <- make_ecoregion(env = rst, nclass = 200)

# Downloading in the European Alps the observations of one plant species
obs.arcto <- get_gbif(
	sp_name = "Arctostaphylos alpinus",
	geo = shp.lonlat,
	grain = 1
)

# Create the range map based on our custom ecoregion at 5 x 5 km resolution
range.arcto <- get_range(
	occ_coord = obs.arcto,
	bioreg = my.eco,
	bioreg_name = "EcoRegion",
	res = 0.05
)

# Plot
countries <- rnaturalearth::ne_countries(
	type = "countries",
	returnclass = "sv"
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
	obs.arcto[, c("decimalLongitude","decimalLatitude")],
	pch = 20,
	col = "#99340470",
	cex = 1
)

###########################################
### Package manuscript plot (Fig 1, arcto)
###########################################

# Root and package
root_dir <- list.files(
	system.file(package = "gbif.range"),
	pattern = "extdata",
	full.names = TRUE
)
if (!dir.exists(file.path(root_dir, "fig_plots"))) {
    dir.create(file.path(root_dir, "fig_plots"))
}
if (!requireNamespace("colorspace", quietly = TRUE)) {
  install.packages("colorspace")
}

# Assign colors to ecoregions
eco.alps <- terra::crop(my.eco, terra::ext(shp.lonlat))
eco.alps2 <- terra::intersect(eco.alps, alps.shp)
col.palette <- colorRampPalette(c("#a6cee3", "#1f78b4",
	"#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
	"#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
	"#ffff99", "#b15928"))
colcol <- col.palette(length(my.eco))
set.seed(7)
eco.alps2$color <- sample(
	paste0(colcol, ""),
	length(eco.alps2),
	replace = FALSE
)

# Extract ecoregions values for points
pt.col <- terra::extract(
	x = eco.alps2,
	y = as.data.frame(obs.arcto[, c("decimalLongitude","decimalLatitude")])
)
pt.plot <- obs.arcto[!is.na(pt.col$color),
						c("decimalLongitude","decimalLatitude")]
pt.col2 <- pt.col[!is.na(pt.col$color), "color"]
pt.col3 <- colorspace::darken(pt.col2, amount=0.6)

# Plot
png(
	paste0(
		root_dir,
		"/fig_plots/fig1_arcto.png"
	),
	width = 100,
	height = 70,
	unit = "cm",
	res = 100,
	pointsize = 110
)
par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 10, cex = 1)
terra::plot(
	eco.alps2,
	col = paste0(eco.alps2$color, "99"),
	border = NA,
	axes = FALSE
)
terra::plot(
	terra::as.polygons(r.arcto),
	border = "black",
	lwd = 7,
	col = "#00000099",
	add = TRUE
)
terra::points(pt.plot, col = pt.col2, pch = 16, cex = 0.7)
terra::points(pt.plot, col = pt.col3, pch = 16, cex = 0.4)
dev.off()

}