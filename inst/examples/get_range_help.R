\dontrun{
###########################################
### Example plot
###########################################

# Load available ecoregions
eco.terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = NULL
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Make range from occurance points
range.tiger <- get_range(
    occ_coord = obs.pt,
    ecoreg = eco.terra,
    ecoreg_name = "ECO_NAME"
)

# Plot
    # Plot political world boundaries
countries <- rnaturalearth::ne_countries(
    type = "countries",
    returnclass = "sv"
)
terra::plot(
    terra::crop(countries, terra::ext(range.tiger$rangeOutput)),
    col = "#bcbddc"
)

    # Plot range 
terra::plot(
    range.tiger$rangeOutput,
    axes = FALSE,
    box = FALSE,
    legend = FALSE,
    col = "chartreuse4",
    add = TRUE
)

    # Plot the occurence points
graphics::points(
    obs.pt[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = "#99340470",
    cex = 1.5
)

###########################################
### Manuscript plot
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
if (!requireNamespace("DescTools", quietly = TRUE)) {
    install.packages("DescTools")
}

# CRS
robin = paste(
    "+proj=robin +lon_0=0 +x_0=0 +y_0=0",
    "+datum=WGS84 +units=m +no_defs +type=crs"
)

# ------------------------------------------------
# Package manuscript plot (Fig 1a, right)
# ------------------------------------------------

# Extent adjust
ext.temp <- terra::ext(range.tiger$rangeOutput)
ext.temp <- c(ext.temp[1]-2, ext.temp[2]+2, ext.temp[3]-2, ext.temp[4]+2)

# Assign colors to ecoregions
eco.local <- terra::crop(eco.terra, ext.temp)
col.palette <- colorRampPalette(
    c("#a6cee3", "#1f78b4", "#b2df8a",
    "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
    "#6a3d9a", "#ffff99", "#b15928")
)
colcol <- col.palette(length(eco.local))
set.seed(3)
eco.local$color <- sample(
    paste0(colcol, ""),
    length(eco.local),
    replace = FALSE
)

# Extract ecoregions values for points
pt.col <- terra::extract(
    eco.local,
    as.data.frame(obs.pt[, c("decimalLongitude","decimalLatitude")])
)
pt.plot = obs.pt[!is.na(pt.col$color),
    c("decimalLongitude","decimalLatitude")]
pt.col2 = pt.col[!is.na(pt.col$color), "color"]
pt.col3 = colorspace::darken(pt.col2, amount=0.6)

# Out points
out.plot <- terra::extract(
    range.tiger$rangeOutput,
    as.data.frame(obs.pt[, c("decimalLongitude","decimalLatitude")])
)
op.na <- is.na(out.plot[,2])
out.plot <- obs.pt[op.na, c("decimalLongitude","decimalLatitude")]

# Plot
png(
    paste0(
        root_dir,
        "/fig_plots/fig1_tiger.png"
    ),
    width = 100,
    height = 70,
    unit = "cm",
    res = 100,
    pointsize = 110
)
par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 10, cex = 0.5)
terra::plot(
    eco.local,
    col = eco.local$color,
    border = NA,
    axes = FALSE
)
terra::plot(
    terra::as.polygons(range.tiger$rangeOutput),
    border = "black",
    lwd = 7,
    col = "#63636399",
    add = TRUE
)
terra::points(pt.plot, col = pt.col2, pch = 16, cex = 1)
terra::points(pt.plot, col = pt.col3, pch = 16, cex = 0.6)
terra::points(out.plot, col = "black", pch = 4, cex = 1, lwd = 10)
dev.off()

# ------------------------------------------------
# Package manuscript plot (Fig 3, world maps)
# ------------------------------------------------

# Load diversity rasters
iucn.path <- paste0(
    system.file(package = "gbif.range"),
    "/extdata/iucn_div_robin.tif"
)
gf.path <- paste0(
    system.file(package = "gbif.range"),
    "/extdata/gf_div_robin.tif"
)
iucn.robin <- terra::rast(iucn.path)
gf.robin <- terra::rast(gf.path)

# Reproject countries
countries.robin <- terra::project(countries, robin)

# Boundary box
bb <- terra::vect(
    rnaturalearth::ne_download(
        type = "wgs84_bounding_box",
        returnclass = "sf",
        category = "physical"
    )
)
bb.robin <- terra::project(bb, robin)

# Plot IUCN map
max.iucn <- round(terra::minmax(iucn.robin)[2])
png(
    paste0(
        root_dir,
        "/fig_plots/fig3_div_iucn.png"
    ),
    width = 100,
    height = 100,
    unit = "cm",
    res = 80,
    pointsize = 80
)
par(mfrow = c(1, 1), lwd = 14, cex = 1.1)
col = colorRampPalette(
    c("#67a9cf", "#f7f7f7", "#ef8a62")
)
terra::plot(
    iucn.robin,
    axes = FALSE,
    legend = FALSE,
    col = col(10),
    smooth = TRUE,
    mar = c(1,1,1,5)
)
terra::plot(countries.robin, add = TRUE, lwd = 5)
terra::plot(bb.robin, add = TRUE, lwd = 8)
par(xpd = NA)
cscl(
    colors = col(10),
    crds = c(-9501111, 9734033, -11800000, -13000000),
    zrng = c(0, max.iucn),
    tickle = -0.3,
    cx = 1.1,
    lablag = -1.3,
    tria = "b",
    at = seq(0, max.iucn, 10),
    horiz = TRUE,
    title = "IUCN richness",
    labs = seq(0, max.iucn, 10),
    titlag = 2
)
par(xpd = FALSE)
dev.off()

# Plot GBIF.RANGE map
max.gf <- round(terra::minmax(gf.robin)[2])
png(
    paste0(
        root_dir,
        "/fig_plots/fig3_div_gf.png"
    ),
    width = 100,
    height = 100,
    unit = "cm",
    res = 80,
    pointsize = 80
)
par(mfrow = c(1, 1), lwd = 14, cex = 1.1)
col = colorRampPalette(
    c("#67a9cf", "#f7f7f7", "#ef8a62")
)
terra::plot(
    gf.robin,
    axes = FALSE,
    legend = FALSE,
    col = col(10),
    smooth = TRUE,
    mar = c(1,1,1,5)
)
terra::plot(countries.robin, add = TRUE, lwd = 5)
terra::plot(bb.robin, add = TRUE, lwd = 8)
par(xpd = NA)
cscl(
    colors = col(10),
    crds = c(-9501111, 9734033, -11800000, -13000000),
    zrng = c(0, max.gf),
    tickle = -0.3,
    cx = 1.1,
    lablag = -1.3,
    tria = "b",
    at = seq(0, max.gf, 10),
    horiz = TRUE,
    title = "gbif.range richness",
    labs = seq(0, max.gf, 10),
    titlag = 2
)
par(xpd = FALSE)
dev.off()

# ------------------------------------------------
# Package manuscript plot (Fig 3, scatter plots)
# ------------------------------------------------

# Load area table
data(area_data)

# Extract data to plot
cor.ras <- rast(list(gf.robin,iucn.robin))
names(cor.ras) <- c("RANGE","IUCN")
set.seed(1)
samp.div <- terra::spatSample(cor.ras,5000,replace=FALSE,na.rm=TRUE)
dat.plot <- list(samp.div,area_data[,-1])
strings <- c("richness","areas")

# Plot diversities relationship
lapply(1:length(dat.plot), function(x) {

    # Extract the data
    xy <- dat.plot[[x]]
    sex <- 2.5
    col <- "#d95f0240"
    add <- ""

    # Log or not depending on the output
    if (strings[x] == "areas") {
        xy <- log(xy + 1)
        sex <- 2.5
        col <- "#d95f0240"
        add <- "(log)"
    }

    # Prep' plotting
    pdf.path <- paste0(
        root_dir,
        "/fig_plots/fig3_",
        strings[x],
        "_gbifrange_cor.pdf"
    )
    pdf(pdf.path, width = 6.3, height = 6)
    par(mfrow = c(1, 1), mar = c(5, 5, 5, 5), lwd = 2)

    # Plot points
    plot(
        xy[, 2], xy[, 1],
        cex.axis = 1.7,
        cex = sex,
        col = col,
        xlim = c(min(xy[, 2]), max(xy[, 2])),
        ylim = c(min(xy[, 1]), max(xy[, 1])),
        pch = 16,
        xlab = paste("IUCN", strings[x], add),
        ylab = paste("gbif.range", strings[x], add),
        cex.lab = 1.7, font.lab = 1
    )

    # Run Linear Models and spearman's correlation
    lm.div <- lm(xy[, 1] ~ xy[, 2])
    ccc.div <- DescTools::CCC(xy[, 2], xy[, 1])$rho.c$est
    cor.div <- cor(xy[, 2],xy[, 1])
    adj.r2 <- summary(lm.div)[[9]]

    # Plot text
    text_cor1 <- bquote("ccc"==.(round(ccc.div,2)))
    text_cor2 <- bquote("r"==.(round(cor.div,2)))
    fig_label(
        text_cor1,
        region = "plot",
        pos = "topleft",
        bty = "n",
        font = 2,
        col = "#121212",
        cex = 2,
        margin = 0.02
    )
    fig_label(
        text_cor2,
        region = "plot",
        pos = "bottomright",
        bty = "n",
        font = 2,
        col = "#6f69c2",
        cex = 2,
        margin = 0.02
    )

    # Plot relationship
    lines(xy[, 2], lm.div$fit, lwd = 7, col = "#7570b3")
    abline(a = 0, b = 1, col = "#252525", lwd = 5, lty = 2)

    # Save
    dev.off()
})

}
