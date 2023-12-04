\dontrun{
# Load available ecoregions
data(ecoregions)

# First download the worldwide observations of Panthera tigris and convert to SpatialPoints
# get occurance points from GBIF
occ  <- get_gbif("Panthera tigris", basis=c("OBSERVATION","HUMAN_OBSERVATION","MACHINE_OBSERVATION"))
# make range from occurance points
range <- get_range("Panthera tigris", occ ,eco.earth ,"ECO_NAME")
# Plot
plot(range, axes = FALSE, box = FALSE, legend=FALSE, col="chartreuse4")
# plot political world boundaries
plot(ne_countries(returnclass = "sf")[1], add=T, col=NA)

# plot(countries,col="#bcbddc")
# plot(ne_countries(returnclass = "sf")[4] , col="antiquewhite")
# TODO FIX THIS
# points(occ[,c("decimalLongitude","decimalLatitude")],pch=20,col="#99340470",cex=1.5)
}
