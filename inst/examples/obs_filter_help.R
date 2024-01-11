
# Load the European Alps extent and a raster of a random resolution
data(geo_dat)
data(exrst)

#read example shapefile
shp.lonlat <- terra::vect("./data/shp_lonlat.shp")

# Downloading in the European Alps the observations of two plant species
obs.arcto = get_gbif("Arctostaphylos alpinus",geo=shp.lonlat)
obs.saxi = get_gbif("Saxifraga cernua",geo=shp.lonlat)
plot(vect(shp.lonlat))
points(obs.arcto[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=1)
points(obs.saxi[,c("decimalLongitude","decimalLatitude")],pch=20,col="#99000d50",cex=1)

# rbind both datasets
both.sp = rbind(obs.arcto,obs.saxi)

# Run function
obs.filt = obs_filter(both.sp,rst)

# Check new points
x11();plot(vect(shp.lonlat))

points(obs.filt[obs.filt$Species%in%"Arctostaphylos alpinus",c("x","y")],
pch=20,col="#238b4550",cex=1)

points(obs.filt[obs.filt$Species%in%"Saxifraga cernua",c("x","y")],
pch=20,col="#99000d50",cex=1)