\donttest{
# Download all worldwide observations of the great panda, with:
# - 100km grain
# - after 1990
# - keeping duplicates
# - adding the name of the person who collected the panda records
obs_am <- get_gbif(
       sp_name = "Ailuropoda melanoleuca",
       grain = 100,
       duplicates = TRUE,
       time_period = c(1990,3000),
       add_infos = c("recordedBy","issue")
)

# Guard against an unavailable remote source
if (gbif_have(obs_am)) {

# Extract borders
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)

# Plot
terra::plot(countries, col = "#bcbddc")
graphics::points(
       obs_am[,c("decimalLongitude","decimalLatitude")],
       pch = 20,
       col = "#238b4550",
       cex = 4
)

}
}
