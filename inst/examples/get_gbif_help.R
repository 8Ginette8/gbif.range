\donttest{
# Download all worldwide observations of Ailuropoda melanoleuca, with:
# - 100km grain
# - after 1990
# - keeping duplicates and by
# - adding the name of the person who collected the panda records
obs.am <- get_gbif(
       sp_name = "Ailuropoda melanoleuca",
       grain = 100 ,
       duplicates = TRUE,
       time_period = c(1990,3000),
       add_infos = c("recordedBy","issue")
)

# Extract borders
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)

# Plot
terra::plot(countries, col = "#bcbddc")
graphics::points(
       obs.am[,c("decimalLongitude","decimalLatitude")],
       pch = 20,
       col = "#238b4550",
       cex = 4
)

}
