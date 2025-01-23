\dontrun{
# Load available ecoregions
eco_terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# First download the worldwide observations of Panthera tigris from GBIF
occ <- get_gbif("Panthera tigris",
				  basis = c("OBSERVATION", 
				          "HUMAN_OBSERVATION", 
				          "MACHINE_OBSERVATION"))

# Make range from occurance points
range <- get_range(occ, eco_terra,"ECO_NAME",clustered_points_outlier = 4)


# Plot

# Plot political world boundaries
countries <- terra::vect(rnaturalearth::ne_countries(returnclass = "sf")[1])
terra::plot(terra::crop(countries,terra::ext(range$range_output)), col = "#bcbddc")

# Plot range 
terra::plot(range$range_output, axes = FALSE, box = FALSE, legend=FALSE,
	col = "chartreuse4", add = TRUE)

# Plot the occurance points
graphics::points(occ[,c("decimalLongitude","decimalLatitude")],pch = 20,col = "#99340470",cex = 1.5)
}
