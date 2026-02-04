\dontrun{
###########################################
### Example plot
###########################################

# Load available ecoregions
eco.terra <- read_bioreg(
    bioreg_name = "eco_terra",
    save_dir = NULL
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Make range from occurance points
range.tiger <- get_range(
    occ_coord = obs.pt,
    bioreg = eco.terra,
    bioreg_name = "ECO_NAME"
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
### Package manuscript plot (Fig 1a, right)
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
    as.data.frame(obs.pt[,c("decimalLongitude","decimalLatitude")])
)
pt.plot = obs.pt[!is.na(pt.col$color),
                    c("decimalLongitude","decimalLatitude")]
pt.col2 = pt.col[!is.na(pt.col$color), "color"]
# install.packages("colorspace")
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

}
