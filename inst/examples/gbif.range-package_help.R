# A self-contained workflow using synthetic environmental data
r1 <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
terra::values(r1) <- rep(seq(0, 1, length.out = 20), each = 20)
r2 <- terra::rast(r1)
terra::values(r2) <- rep(seq(0, 1, length.out = 20), times = 20)
env <- c(r1, r2)

# Create custom ecoregions and a small occurrence dataset
eco <- make_ecoreg(env = env, nclass = 4)
occ <- data.frame(
  decimalLongitude = c(0.5, 1.2, 2.4, 3.6, 2.8, 6.0, 7.2, 7.4, 7.6, 8.8),
  decimalLatitude = c(1.0, 0.1, 2.3, 2.5, 2.7, 5.0, 7.1, 7.3, 6.5, 7.7)
)

# Build a range map
range_obj <- get_range(
  occ_coord = occ,
  ecoreg = eco,
  verbose = FALSE)

# Plot the range map
terra::plot(range_obj$rangeOutput, col=3, main = "Range Map")
# Plot the occurrence points
points(occ, pch="x")

\dontrun{
# Typical online workflow with GBIF data
obs <- get_gbif("Panthera tigris", grain = 25)
# obs <- obs[obs$taxonRank=="SPECIES",]
# status <- get_status("Panthera tigris")
eco_terra <- read_ecoreg("eco_terra")
tiger_range <- get_range(
  occ_coord = obs,
  ecoreg = eco_terra,
  ecoreg_name = "ECO_NAME"
)
plot(tiger_range$rangeOutput, col=3, main = paste("Range:", obs$scientificName[1]))
points(obs$decimalLongitude, obs$decimalLatitude,  pch="x", col=rgb(1, 0, 1, 0.2))
}
