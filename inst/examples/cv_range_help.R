\dontrun{
###########################################
### Example plot
###########################################

# Load available ecoregions
eco.terra <- read_bioreg(
    bioreg_name = "eco_terra",
    save_dir = NULL,
    format = "sf"
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Make range from occurance points
range.tiger <- get_range(
    occ_coord = obs.pt,
    bioreg = eco.terra,
    bioreg_name = "ECO_NAME",
    format = "sf"
)

pt.test <- cv_range(
    range_object = range.tiger,
    cv = 'block-cv',
    nfolds = 5,
    nblocks = 2
);pt.test


###########################################
### Package manuscript plot (Fig 2a-b)
###########################################

# Root and package
root_dir <- list.files(
    system.file(package = "gbif.range"),
    pattern = "extdata",
    full.names = TRUE
)
if (!dir.exists(file.path(root_dir, "fig_plots"))) {
    dir.create(file.path(root_dir, "fig_plots"))
}
if (!requireNamespace("colorspace", quietly = TRUE)) {
  install.packages("colorspace")
}

###########
##### Plant #########
###########

# Preliminary
spdf.world <- terra::vect(
    rnaturalearth::ne_countries(scale=10,returnclass="sf")
)
shp.lonlat <- terra::vect(
    paste0(
        system.file(package = "gbif.range"),
        "/extdata/shp_lonlat.shp"
    )
)
obs.arcto <- get_gbif("Arctostaphylos alpinus", geo = shp.lonlat, grain = 1)
rst <- terra::rast(
    paste0(
        system.file(package = "gbif.range"),
        "/extdata/rst.tif"
    )
)
my.eco <- make_ecoregion(rst, 200)
range.arcto <- get_range(
    occ_coord = obs.arcto,
    bioreg = my.eco,
    bioreg_name = "EcoRegion",
    res = 0.05
)
ext.temp <- terra::ext(range.arcto$rangeOutput)
ext.temp <- c(ext.temp[1]-0.6, ext.temp[2]+0.05,
                ext.temp[3]-0.05, ext.temp[4]+0.05)

# Create pseudo-absences
    # Firt remove observations considered as outliers in get_range
xy.df <- range.arcto$init.args$occ_coord
r.ext <- terra::ext(range.arcto$rangeOutput)
Xrm.cond <- xy.df$decimalLongitude >= r.ext[1] &
                    xy.df$decimalLongitude <= r.ext[2]
Yrm.cond <-  xy.df$decimalLatitude >= r.ext[3] &
                    xy.df$decimalLatitude <= r.ext[4]
xy.df <- xy.df[Xrm.cond & Yrm.cond,]

    # Samples n regular background points over the original range extent
x.interv <- (r.ext[2] - r.ext[1]) / (sqrt(1e4)-1)
y.interv <- (r.ext[4] - r.ext[3]) / (sqrt(1e4)-1)
lx <- seq(r.ext[1], r.ext[2], x.interv)
ly <- seq(r.ext[3], r.ext[4], y.interv)
bp.xy <- expand.grid(decimalLongitude = lx, decimalLatitude = ly)

    # Combine observations with background
obs.xy <- xy.df[, c("decimalLongitude","decimalLatitude")]
all.xy <- rbind(obs.xy, bp.xy)
all.xy$Pres <- 0
all.xy[1:nrow(obs.xy), "Pres"] <- 1

    # Run block-cv
xy.pres <- all.xy$Pres
cv.strat <- make_blocks(
    nfolds = 5,
    df = all.xy[, c("decimalLongitude","decimalLatitude")],
    nblocks = 5*2,
    pres = xy.pres
)
all.xy$bcv <- cv.strat
all.xy[all.xy$bcv%in%1, "col"] <- "#e41a1c"
all.xy[all.xy$bcv%in%2, "col"] <- "#377eb8"
all.xy[all.xy$bcv%in%3, "col"] <- "#4daf4a"
all.xy[all.xy$bcv%in%4, "col"] <- "#984ea3"
all.xy[all.xy$bcv%in%5, "col"] <- "#ff7f00"
all.xy[all.xy$Pres%in%1, "col"] <- colorspace::darken(
    all.xy[all.xy$Pres%in%1, "col"],
    amount = 0.5
)

# Plot
    # Evaluate
ar.test <- cv_range(
    range_object = range.arcto,
    cv = 'block-cv',
    nfolds = 5,
    nblocks = 2
)

    # First extract the world at the arcto extent and divide
    # presences/absences/outliers
world.local <- terra::crop(spdf.world, ext.temp)
world.local.ar <- terra::aggregate(world.local)
pres <- all.xy[all.xy$Pres%in%1, ]
abs <- all.xy[all.xy$Pres%in%0, ]
id.in <- terra::extract(
    range.arcto$rangeOutput,
    as.data.frame(pres[,1:2])
)
pres <- pres[!is.na(id.in$layer), ]

    # Continue plotting
png(
    paste0(
        root_dir,
        "/fig_plots/fig2_arcto_cv.png"
    ),
    width = 100,
    height = 70,
    unit = "cm",
    res = 100,
    pointsize = 110
)
par(
    mfrow = c(1,1),
    mar = c(5,5,5,20),
    lwd = 1,
    cex = 0.5
)
terra::plot(world.local.ar, col = "#bcd1bc", axes = FALSE, lwd = 2)
terra::points(abs, col = paste0(abs$col, "50"), pch = 16, cex = 0.7)
terra::plot(
    terra::as.polygons(range.arcto$rangeOutput),
    border = "black",
    lwd = 6,
    col = "#63636370",
    add = TRUE
)
terra::points(pres, col = paste0(pres$col,"90"), pch = 16, cex = 1.5)
text(
    6.6,43.2,
    paste("Mean TSS =", round(tail(ar.test[,"TSS"],1),2)),
    cex = 1.5,
    font = 2
)
text(
    7.2,42.8,
    paste("Mean Precision =", round(tail(ar.test[,"Precision"],1),2)),
    cex = 1.5,
    font = 2
)
dev.off()

###########
##### Tiger #########
###########

# Preliminary
ext.temp <- terra::ext(range.tiger$rangeOutput)
ext.temp <- c(ext.temp[1]-0.2, ext.temp[2]+0.2,
                ext.temp[3]-0.2, ext.temp[4]+0.2)

# Create pseudo-absences
    # Firt remove observations considered as outliers in get_range
xy.df <- range.tiger$init.args$occ_coord
r.ext <- terra::ext(range.tiger$rangeOutput)
Xrm.cond <- xy.df$decimalLongitude >= r.ext[1] &
                xy.df$decimalLongitude <= r.ext[2]
Yrm.cond <-  xy.df$decimalLatitude >= r.ext[3] &
                xy.df$decimalLatitude <= r.ext[4]
xy.df <- xy.df[Xrm.cond & Yrm.cond,]

    # Samples n regular background points over the original range extent
x.interv <- (r.ext[2] - r.ext[1]) / (sqrt(1e4) - 1)
y.interv <- (r.ext[4] - r.ext[3]) / (sqrt(1e4) - 1)
lx <- seq(r.ext[1], r.ext[2], x.interv)
ly <- seq(r.ext[3], r.ext[4], y.interv)
bp.xy <- expand.grid(decimalLongitude = lx, decimalLatitude = ly)

    # Combine observations with background
obs.xy <- xy.df[,c("decimalLongitude", "decimalLatitude")]
all.xy <- rbind(obs.xy, bp.xy)
all.xy$Pres <- 0
all.xy[1:nrow(obs.xy), "Pres"] <- 1

    # Run block-cv
xy.pres <- all.xy$Pres
cv.strat <- make_blocks(
    nfolds = 5,
    df = all.xy[, c("decimalLongitude","decimalLatitude")],
    nblocks = 5*2,
    pres = xy.pres
)
all.xy$bcv <- cv.strat
all.xy[all.xy$bcv%in%1,"col"] <- "#e41a1c"
all.xy[all.xy$bcv%in%2,"col"] <- "#377eb8"
all.xy[all.xy$bcv%in%3,"col"] <- "#4daf4a"
all.xy[all.xy$bcv%in%4,"col"] <- "#984ea3"
all.xy[all.xy$bcv%in%5,"col"] <- "#ff7f00"
all.xy[all.xy$Pres%in%1,"col"] <- colorspace::darken(
    all.xy[all.xy$Pres%in%1, "col"],
    amount = 0.5
)

# Plot
    # First extract the world at the tiger extent and divide
    # presences/absences/outliers
world.local <- terra::crop(spdf.world, ext.temp)
world.local.ti <- terra::aggregate(world.local)
pres <- all.xy[all.xy$Pres%in%1, ]
abs <- all.xy[all.xy$Pres%in%0, ]
id.in <- terra::extract(
    range.tiger$rangeOutput,
    as.data.frame(pres[,1:2])
)
pres <- pres[!is.na(id.in$layer),]

    # Continue plotting
png(
    paste0(
        root_dir,
        "/fig_plots/fig2_tiger_cv.png"
    ),
    width = 100,
    height = 70,
    unit = "cm",
    res = 100,
    pointsize = 110
)
par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 1, cex = 0.5)
terra::plot(world.local.ti, col = "#bcd1bc", axes = FALSE, lwd = 2)
terra::points(abs, col = paste0(abs$col,"50"), pch = 16, cex = 0.6)
terra::plot(
    terra::as.polygons(range.tiger$rangeOutput),
    border = "black",
    lwd = 6,
    col = "#63636370",
    add = TRUE
)
terra::points(pres, col = paste0(pres$col,"80"), pch = 16, cex = 1.6)
text(86,49,
    paste("Mean TSS =", round(tail(pt.test[,"TSS"],1),2)),
    cex = 1.5,
    font = 2
)
text(90.4,45.4,
    paste("Mean Precision =", round(tail(pt.test[,"Precision"],1),2)),
    cex = 1.5, font = 2
)
dev.off()

}
