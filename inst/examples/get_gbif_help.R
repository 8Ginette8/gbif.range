# Downloading worldwide the observations of Panthera tigris
obs.pt <- get_gbif("Panthera tigris",basis=c("OBSERVATION","HUMAN_OBSERVATION","MACHINE_OBSERVATION"))
countries <- terra::vect(rnaturalearth::ne_countries(type = "countries",returnclass = "sf"))
terra::plot(countries,col="#bcbddc")
points(obs.pt[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=4)

\dontrun{
# Downloading worldwide the observations of Ailuropoda melanoleuca (with a 100km grain, after 1990
# and by keeping duplicates and by adding the name of the person who collected the panda records)
obs.am <- get_gbif("Ailuropoda melanoleuca", grain = 100000 , duplicates = TRUE,
   time_period = c(1990,3000), add_infos = c("recordedBy","issue"))
terra::plot(countries,col="#bcbddc")
graphics::points(obs.am[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=4)

# Downloading worlwide the observations of Phascolarctos cinereus (with a 1km grain, after 1980,
# and keeping raster centroids)
obs.pc <- get_gbif("Phascolarctos cinereus", grain = 1000,
   time_period = c(1990,3000), centroids = TRUE)
}