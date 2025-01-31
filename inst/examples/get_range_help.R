\dontrun{
# Load available ecoregions
eco_terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# First download the worldwide observations of Panthera tigris from GBIF
occ <- get_gbif(sp_name = "Panthera tigris",
                time_period = c(2000, 3000),
                basis = c("OBSERVATION",
                          "HUMAN_OBSERVATION",
                          "MACHINE_OBSERVATION",
                          "OCCURRENCE"))

# Make range from occurance points
range <- get_range(occ, eco_terra, "ECO_NAME", clustered_points_outlier = 4)


# Plot

# Plot political world boundaries
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(terra::crop(countries, terra::ext(range$range_output)), col = "#bcbddc")

# Plot range 
terra::plot(range$range_output, axes = FALSE, box = FALSE, legend = FALSE,
	col = "chartreuse4", add = TRUE)

# Plot the occurance points
graphics::points(occ[, c("decimalLongitude","decimalLatitude")],
    pch = 20,col = "#99340470",cex = 1.5)
}
