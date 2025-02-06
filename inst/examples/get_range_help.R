\dontrun{
# Load available ecoregions
eco.terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pt <- get_gbif(sp_name = "Panthera tigris",
                   basis = c("OBSERVATION",
                          "HUMAN_OBSERVATION",
                          "MACHINE_OBSERVATION",
                          "OCCURRENCE"))

# Make range from occurance points
range.tiger <- get_range(occ_coord = obs.pt,
                         bioreg = eco.terra,
                         bioreg_name = "ECO_NAME",
                         degrees_outlier = 6,
                         clustered_points_outlier = 4)

# Plot

# Plot political world boundaries
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(terra::crop(countries, terra::ext(range.tiger$rangeOutput)), col = "#bcbddc")

# Plot range 
terra::plot(range.tiger$rangeOutput, axes = FALSE, box = FALSE, legend = FALSE,
	col = "chartreuse4", add = TRUE)

# Plot the occurance points
graphics::points(obs.pt[, c("decimalLongitude","decimalLatitude")],
    pch = 20,col = "#99340470",cex = 1.5)
}
