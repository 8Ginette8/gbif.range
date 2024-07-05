#TODO transfer to data()
rst.path <- paste0(system.file(package = "gbif.range"),"/extdata/rst_enl.tif")
rst <- terra::rast(rst.path)
shp.path <- paste0(system.file(package = "gbif.range"),"/extdata/shp_lonlat.shp")
shp.lonlat <- terra::vect(shp.path)
rst <- terra::crop(rst, shp.lonlat)
#plot(crp)

# Apply the function by infering 50 classes of environments
my.eco <- make_ecoregion(rst,50)
terra::plot(my.eco)

# Downloading in the European Alps the observations of one plant species
obs.arcto <- get_gbif("Arctostaphylos alpinus",geo=shp.lonlat)

# Create the range map based on our custom ecoregion
range.arcto <- get_range(obs.arcto, my.eco, "EcoRegion", res=20)

# Plot
terra::plot(shp.lonlat, col="grey")
terra::plot(range.arcto,add=TRUE,col=rgb(0.2,1,0.2,0.5,1))
graphics::points(obs.arcto[,c("decimalLongitude","decimalLatitude")],pch=20,col="orange1",cex=1)