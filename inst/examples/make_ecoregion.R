# Load the European Alps extent and a raster of a random resolution
data(geo_dat)
data(exrst)

# Apply the function by infering 50 classes of environments
my.eco = make_ecoregion(rst,50)
plot(my.eco)

# Downloading in the European Alps the observations of one plant species
obs.arcto = get_gbif("Arctostaphylos alpinus",geo=shp.lonlat)

# Create the range map based on our custom ecoregion
range.arcto = get_range("Arctostaphylos alpinus",obs.arcto,my.eco,"EcoRegion",res=20)

# Plot
plot(vect(shp.lonlat))
plot(range.arcto,add=TRUE,col="darkgreen")
points(obs.arcto[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=1)