\donttest{
# Downloading worldwide all the observations of the great panda
obs_am <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Guard against an unavailable remote source
if (gbif_have(obs_am)) {

# Create a vector of folds (n = 5) spatially blocked (n = 10)
block_am <- make_blocks(
    nfolds = 5,
    df = obs_am[, c("decimalLatitude","decimalLongitude")],
    nblocks = 5
)

# Plot one colour per fold
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
countries_focus <- terra::crop(
    countries,
    terra::ext(73.0, 135.0, 18.0, 54.0)
)
terra::plot(countries_focus, col = "#bcbddc")
graphics::points(
    obs_am[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = block_am,
    cex = 1
)

}
}
