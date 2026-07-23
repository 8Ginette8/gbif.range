\donttest{
# Downloading worldwide all the observations of great panda 
obs.pd <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Create a vector of folds (n = 5) spatially blocked (n = 10)
block.pd <- make_blocks(
    nfolds = 5,
    df = obs.pd[, c("decimalLatitude","decimalLongitude")],
    nblocks = 5
)

# Plot one colour per fold
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
countries.focus <- terra::crop(
    countries,
    terra::ext(73.0, 135.0, 18.0, 54.0)
)
terra::plot(countries.focus, col = "#bcbddc")
graphics::points(
    obs.pd[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = block.pd,
    cex = 1
)

}
