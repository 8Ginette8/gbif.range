\dontrun{
# Load available ecoregions
eco_terra = read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# First download the worldwide observations of Panthera tigris from GBIF
occ <- get_gbif("Panthera tigris",
				  basis=c("OBSERVATION","HUMAN_OBSERVATION","MACHINE_OBSERVATION"))

# Make range from occurance points
range <- get_range("Panthera tigris", occ ,eco_terra ,"ECO_NAME")

# Plot
terra::plot(range, axes = FALSE, box = FALSE, legend=FALSE, col="chartreuse4")

# Plot political world boundaries
terra::plot(rnaturalearth::ne_countries(returnclass = "sf")[1], add=T, col=NA)
graphics::points(occ[,c("decimalLongitude","decimalLatitude")],pch=20,col="#99340470",cex=1.5)

}
