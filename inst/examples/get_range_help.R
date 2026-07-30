\donttest{
# Load available ecoregions
eco_terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir()
)

# First download the whole available observations of the panda from GBIF
obs_am <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Guard against an unavailable remote source
if (gbif_have(eco_terra, obs_am)) {

# Build a range map from occurrence points
range_panda <- get_range(
    occ_coord = obs_am,
    ecoreg = eco_terra,
    ecoreg_name = "ECO_NAME",
    format = "SpatRaster"
)

# Plot
    # Plot political world boundaries
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
terra::plot(
    terra::crop(countries, terra::ext(73.0, 135.0, 18.0, 54.0)),
    col = "#bcbddc"
)

    # Plot range
terra::plot(
    range_panda$rangeOutput,
    axes = FALSE,
    box = FALSE,
    legend = FALSE,
    col = "chartreuse4",
    main = paste("Range:", obs_am$scientificName[1]),
    add = TRUE
)

    # Plot the occurrence points
graphics::points(
    obs_am[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = "#99340470",
    cex = 1.5
)

}
}
