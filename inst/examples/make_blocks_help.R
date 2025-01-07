\dontrun{
# Downloading worldwide the observations of Panthera tigris
obs.pt <- get_gbif("Panthera tigris",
                   basis = c("OBSERVATION",
                           "HUMAN_OBSERVATION",
                           "MACHINE_OBSERVATION"))

# Create a vector of folds (n = 5) spatially blocked (n = 10)
block.pt <- make_blocks(nfolds = 5,
        df = obs.pt[,c("decimalLatitude","decimalLongitude")],
        nblocks = 10)

# Plot one colour per fold
countries <- terra::vect(rnaturalearth::ne_countries(type = "countries", returnclass = "sf"))
countries.focus <- terra::crop(countries, terra::ext(60,100,0,40))
terra::plot(countries.focus, col = "#bcbddc")
graphics::points(obs.pt[,c("decimalLongitude","decimalLatitude")],
        pch = 20, col = block.pt, cex = 1)
}