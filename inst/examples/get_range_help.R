\donttest{
# Load available ecoregions
eco.terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir()
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pd <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Build a range map from occurrence points
range.panda <- get_range(
    occ_coord = obs.pd,
    ecoreg = eco.terra,
    clust_pts_outlier = 6,
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
    range.panda$rangeOutput,
    axes = FALSE,
    box = FALSE,
    legend = FALSE,
    col = "chartreuse4",
    main = paste("Range:", obs.pd$scientificName[1]),
    add = TRUE
)

    # Plot the occurrence points
graphics::points(
    obs.pd[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = "#99340470",
    cex = 1.5
)

}
