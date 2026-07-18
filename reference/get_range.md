# Create Species Range Maps from Occurrences and Ecoregions

Estimate ecologically informed species ranges from occurrence data and
an ecoregion layer. The workflow combines outlier filtering, clustering,
convex hull construction, and intersection with occupied ecoregions.

## Usage

``` r
get_range(
  occ_coord = NULL,
  ecoreg = NULL,
  ecoreg_name = NULL,
  degrees_outlier = 5,
  clust_pts_outlier = 4,
  buff_width_point = 4,
  buff_incrmt_pts_line = 0.5,
  buff_width_polygon = 4,
  dir_temp = tempdir(),
  format = c("SpatVector", "sf", "SpatRaster"),
  res = 0.1,
  verbose = TRUE
)
```

## Arguments

- occ_coord:

  `getGBIF` object returned by
  [`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
  or a `data.frame` containing the columns `decimalLongitude` and
  `decimalLatitude`.

- ecoreg:

  Spatial ecoregion layer in WGS84. Accepted classes are
  `SpatialPolygonsDataFrame`, `SpatVector`, and `sf`. This can be a
  downloaded layer from
  [`read_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)
  or a custom layer created by
  [`make_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md).

- ecoreg_name:

  Character. String naming the field in `ecoreg` that defines ecoregion
  categories. For `eco_terra`, for example, `"ECO_NAME"` is the most
  detailed level. When `ecoreg` comes from
  [`make_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md),
  `"EcoRegion"` is used automatically.

- degrees_outlier:

  Numeric. Distance threshold in degrees. Points whose
  `clust_pts_outlier`-th nearest neighbour lies beyond this distance are
  treated as outliers. Default is `5`.

- clust_pts_outlier:

  Numeric. k-nearest-neighbour order used for outlier detection. Default
  is `4`.

- buff_width_point:

  Numeric. Buffer width in degrees for isolated points.

- buff_incrmt_pts_line:

  Numeric. Increment in buffer width for linear clusters.

- buff_width_polygon:

  Numeric. Buffer width in degrees applied to convex hull polygons.

- dir_temp:

  Character. String giving the directory used for temporary convex-hull
  files. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- format:

  Character. Output format for the range map. One of `"SpatVector"`
  (default), `"sf"`, or `"SpatRaster"`. `"SpatRaster"` rasterizes the
  range at the resolution set by `res`.

- res:

  Numeric. Output resolution in degrees when `format = "SpatRaster"`.
  Default is `0.1` (about 11.1 km at the equator). The achievable
  resolution is constrained by the spatial precision of `ecoreg`.

- verbose:

  Logical. Should progress messages be printed?

## Value

An object of class `getRange` with two fields: `init.args`, containing
the arguments and data used to build the map, and `rangeOutput`,
containing the resulting `SpatVector`, `sf`, or `SpatRaster` object.

## Details

The function follows four main steps.

First, occurrence points are filtered for spatial outliers using
k-nearest-neighbour distances and are assigned to ecoregions.

Second, points within occupied ecoregions are grouped into clusters
using a combination of Gaussian mixture modeling and k-means clustering.

Third, each cluster is converted into a polygon. Single points receive
circular buffers, collinear clusters receive line-based buffers, and
other clusters receive buffered convex hulls through `conv_function()`.

Fourth, cluster polygons are intersected with their parent ecoregions
and merged into a final range layer.

Download-ready ecoregion datasets include `eco_terra` for terrestrial
species, `eco_fresh` for freshwater species, and `eco_marine` or
`eco_hd_marine` for marine species.

## References

Hagen, O., Vaterlaus, L., Albouy, C., Brown, A., Leugger, F., Onstein,
R. E., Novaes de Santana, C., Scotese, C. R., & Pellissier, L. (2019).
Mountain building, climate cooling and the richness of cold-adapted
plants in the Northern Hemisphere. Journal of Biogeography.
[doi:10.1111/jbi.13653](https://doi.org/10.1111/jbi.13653)

Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J.
S., ... & Pellissier, L. (2022). An integrated high resolution mapping
shows congruent biodiversity patterns of Fagales and Pinales. New
Phytologist, 235(2), 759-772.
[doi:10.1111/nph.18158](https://doi.org/10.1111/nph.18158)

Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand,
H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H.,
Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R.
(2001). Terrestrial ecoregions of the world: a new map of life on Earth.
BioScience, 51(11), 933-938.
[doi:10.1641/0006-3568(2001)051\[0933:TEOTWA\]2.0.CO;2](https://doi.org/10.1641/0006-3568%282001%29051%5B0933%3ATEOTWA%5D2.0.CO%3B2)

The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
Units. GIS layers developed by The Nature Conservancy with multiple
partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986.
Cambridge (UK): The Nature Conservancy.

Spalding, M. D., Fox, H. E., Allen, G. R., Davidson, N., Ferdana, Z. A.,
Finlayson, M., Halpern, B. S., Jorge, M. A., Lombana, A., Lourie, S. A.,
Martin, K. D., McManus, E., Molnar, J., Recchia, C. A., Robertson, J.
(2007). Marine Ecoregions of the World: A Bioregionalization of Coastal
and Shelf Areas. BioScience, 57(7), 573-583.
[doi:10.1641/B570707](https://doi.org/10.1641/B570707)

Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
Pelagic provinces of the world: a biogeographic classification of the
world's surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
[doi:10.1016/j.ocecoaman.2011.12.016](https://doi.org/10.1016/j.ocecoaman.2011.12.016)

The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
of the World. GIS layers developed by The Nature Conservancy with
multiple partners, combined from Spalding et al. (2007) and Spalding et
al. (2012). Cambridge (UK): The Nature Conservancy.

Abell, R., Thieme, M. L., Revenga, C., Bryer, M., Kottelat, M.,
Bogutskaya, N., Coad, B., Mandrak, N., Contreras Balderas, S., Bussing,
W., Stiassny, M. L. J., Skelton, P., Allen, G. R., Unmack, P., Naseka,
A., Ng, R., Sindorf, N., Robertson, J., Armijo, E., Higgins, J. V.,
Heibel, T. J., Wikramanayake, E., Olson, D., Lopez, H. L., Reis, R. E.,
Lundberg, J. G., Sabaj Perez, M. H., Petry, P. (2008). Freshwater
Ecoregions of the World: A New Map of Biogeographic Units for Freshwater
Biodiversity Conservation. BioScience, 58(5), 403-414.
[doi:10.1641/B580507](https://doi.org/10.1641/B580507)

Hijmans, R. J. (2022). terra: Spatial Data Analysis. R package version
1.6-7. <https://cran.r-project.org/package=terra>

## See also

[`read_ecoreg`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)()
and
[`make_ecoreg`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md)()
to prepare the ecoregion layer used here;
[`cv_range`](https://8ginette8.github.io/gbif.range/reference/cv_range.md)()
and
[`evaluate_range`](https://8ginette8.github.io/gbif.range/reference/evaluate_range.md)()
to evaluate the resulting range map.

## Examples

``` r
# \donttest{
###########################################
### Example plot
###########################################

# Load available ecoregions
eco.terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir()
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pt <- get_gbif(sp_name = "Panthera tigris")
#> |--------------------------------------------|
#> | Total number (all records)    :       8007 |
#> | Kept records                  :       5457 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE
#> 
#> ...GBIF records of Panthera tigris: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> -----------------------------------------------
#>                     step removed remaining
#>          Grain filtering     117      5340
#>       Duplicated records    2587      2753
#>          Absence records       0      2753
#>          Basis selection      82      2671
#>  Establishment selection       0      2671
#>               Time frame       0      2671
#>        Identical records       0      2671
#>         Raster centroids       0      2671
#> 
#> Initial records         : 5457
#> Total removed           : 2786
#> Final records (XY)      : 2671
#> -----------------------------------------------
#> Final records (no XY)   : 0

# Build a range map from occurrence points
range.tiger <- get_range(
    occ_coord = obs.pt,
    ecoreg = eco.terra,
    ecoreg_name = "ECO_NAME",
    format = "SpatRaster"
)
#> ## Start of computation for species:  Panthera tigris  ### 
#> 27 outlier's from 2666 | proportion from total points: 1%
#> ecoregion 1  of  44 :  Amur Meadow Steppe 
#> ecoregion 2  of  44 :  Brahmaputra Valley Semi-Evergreen Forests 
#> ecoregion 3  of  44 :  Central Deccan Plateau Dry Deciduous Forests 
#> ecoregion 4  of  44 :  Central Indochina Dry Forests 
#> ecoregion 5  of  44 :  Chao Phraya Freshwater Swamp Forests 
#> ecoregion 6  of  44 :  Chhota-Nagpur Dry Deciduous Forests 
#> ecoregion 7  of  44 :  Deccan Thorn Scrub Forests 
#> ecoregion 8  of  44 :  Eastern Highlands Moist Deciduous Forests 
#> ecoregion 9  of  44 :  Eastern Himalayan Alpine Shrub And Meadows 
#> ecoregion 10  of  44 :  Eastern Himalayan Broadleaf Forests 
#> ecoregion 11  of  44 :  Eastern Himalayan Subalpine Conifer Forests 
#> ecoregion 12  of  44 :  Eastern Java-Bali Rain Forests 
#> ecoregion 13  of  44 :  Himalayan Subtropical Broadleaf Forests 
#> ecoregion 14  of  44 :  Himalayan Subtropical Pine Forests 
#> ecoregion 15  of  44 :  Kayah-Karen Montane Rain Forests 
#> ecoregion 16  of  44 :  Khathiar-Gir Dry Deciduous Forests 
#> ecoregion 17  of  44 :  Lower Gangetic Plains Moist Deciduous Forests 
#> ecoregion 18  of  44 :  Luang Prabang Montane Rain Forests 
#> ecoregion 19  of  44 :  Manchurian Mixed Forests 
#> ecoregion 20  of  44 :  Meghalaya Subtropical Forests 
#> ecoregion 21  of  44 :  Mizoram-Manipur-Kachin Rain Forests 
#> ecoregion 22  of  44 :  Narmada Valley Dry Deciduous Forests 
#> ecoregion 23  of  44 :  North Western Ghats Moist Deciduous Forests 
#> ecoregion 24  of  44 :  North Western Ghats Montane Rain Forests 
#> ecoregion 25  of  44 :  Northwestern Thorn Scrub Forests 
#> ecoregion 26  of  44 :  Okhotsk-Manchurian Taiga 
#> ecoregion 27  of  44 :  Orissa Semi-Evergreen Forests 
#> ecoregion 28  of  44 :  Peninsular Malaysian Montane Rain Forests 
#> ecoregion 29  of  44 :  Peninsular Malaysian Rain Forests 
#> ecoregion 30  of  44 :  South Deccan Plateau Dry Deciduous Forests 
#> ecoregion 31  of  44 :  South Western Ghats Moist Deciduous Forests 
#> ecoregion 32  of  44 :  South Western Ghats Montane Rain Forests 
#> ecoregion 33  of  44 :  Southeastern Indochina Dry Evergreen Forests 
#> ecoregion 34  of  44 :  Suiphun-Khanka Meadows And Forest Meadows 
#> ecoregion 35  of  44 :  Sumatran Freshwater Swamp Forests 
#> ecoregion 36  of  44 :  Sumatran Lowland Rain Forests 
#> ecoregion 37  of  44 :  Sumatran Montane Rain Forests 
#> ecoregion 38  of  44 :  Sumatran Peat Swamp Forests 
#> ecoregion 39  of  44 :  Sunda Shelf Mangroves 
#> ecoregion 40  of  44 :  Sundarbans Mangroves 
#> ecoregion 41  of  44 :  Tenasserim-South Thailand Semi-Evergreen Rain Forests 
#> ecoregion 42  of  44 :  Terai-Duar Savanna And Grasslands 
#> ecoregion 43  of  44 :  Upper Gangetic Plains Moist Deciduous Forests 
#> ecoregion 44  of  44 :  Ussuri Broadleaf And Mixed Forests 
#> ## End of computation for species:  Panthera tigris  ### 

# Plot
    # Plot political world boundaries
countries <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
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

    # Plot the occurrence points
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
root_dir <- tempdir()
if (!dir.exists(file.path(root_dir, "fig_plots"))) {
    dir.create(file.path(root_dir, "fig_plots"))
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
col.palette <- grDevices::colorRampPalette(
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
pt.col3 = grDevices::adjustcolor(pt.col2, red.f = 0.6, green.f = 0.6, blue.f = 0.6)

# Out points
out.plot <- terra::extract(
    range.tiger$rangeOutput,
    as.data.frame(obs.pt[, c("decimalLongitude","decimalLatitude")])
)
op.na <- is.na(out.plot[,2])
out.plot <- obs.pt[op.na, c("decimalLongitude","decimalLatitude")]

# Plot
grDevices::png(
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
oldpar <- graphics::par(mfrow = c(1, 1), mar = c(5, 5, 5, 20), lwd = 10, cex = 0.5)
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
graphics::points(pt.plot, col = pt.col2, pch = 16, cex = 1)
graphics::points(pt.plot, col = pt.col3, pch = 16, cex = 0.6)
graphics::points(out.plot, col = "black", pch = 4, cex = 1, lwd = 10)
grDevices::dev.off()
#> agg_record_1ce32c50c774 
#>                       2 
graphics::par(oldpar)

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
bb <- terra::as.polygons(
  terra::ext(-180, 180, -90, 90),
  crs = "EPSG:4326"
)
bb <- terra::densify(bb, interval = 1)
bb.robin <- terra::project(bb, robin)

# Plot IUCN map
max.iucn <- round(terra::minmax(iucn.robin)[2])
grDevices::png(
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
oldpar <- graphics::par(mfrow = c(1, 1), lwd = 14, cex = 1.1)
col = grDevices::colorRampPalette(
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
graphics::par(xpd = NA)
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
grDevices::dev.off()
#> agg_record_1ce32c50c774 
#>                       2 
graphics::par(oldpar)

# Plot GBIF.RANGE map
max.gf <- round(terra::minmax(gf.robin)[2])
grDevices::png(
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
oldpar <- graphics::par(mfrow = c(1, 1), lwd = 14, cex = 1.1)
col = grDevices::colorRampPalette(
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
graphics::par(xpd = NA)
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
grDevices::dev.off()
#> agg_record_1ce32c50c774 
#>                       2 
graphics::par(oldpar)

# ------------------------------------------------
# Package manuscript plot (Fig 3, scatter plots)
# ------------------------------------------------

# Load area table
data(area_data)

# Extract data to plot
cor.ras <- terra::rast(list(gf.robin,iucn.robin))
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
    grDevices::pdf(pdf.path, width = 6.3, height = 6)
    oldpar <- graphics::par(mfrow = c(1, 1), mar = c(5, 5, 5, 5), lwd = 2)

    # Plot points
    graphics::plot(
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
    lm.div <- stats::lm(xy[, 1] ~ xy[, 2])

    # Lin's concordance correlation coefficient (population covariance/variance)
    ccc_lin <- function(x, y) {
      mx <- mean(x); my <- mean(y)
      vx <- mean((x - mx)^2)
      vy <- mean((y - my)^2)
      sxy <- mean((x - mx) * (y - my))
      2 * sxy / (vx + vy + (mx - my)^2)
    }
    ccc.div <- ccc_lin(xy[, 2], xy[, 1])
    cor.div <- stats::cor(xy[, 2],xy[, 1])
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
    graphics::lines(xy[, 2], lm.div$fit, lwd = 7, col = "#7570b3")
    graphics::abline(a = 0, b = 1, col = "#252525", lwd = 5, lty = 2)

    # Save
    grDevices::dev.off()
    graphics::par(oldpar)
})
#> [[1]]
#> [[1]]$mfrow
#> [1] 1 1
#> 
#> [[1]]$mar
#> [1] 5.1 4.1 4.1 2.1
#> 
#> [[1]]$lwd
#> [1] 1
#> 
#> 
#> [[2]]
#> [[2]]$mfrow
#> [1] 1 1
#> 
#> [[2]]$mar
#> [1] 5.1 4.1 4.1 2.1
#> 
#> [[2]]$lwd
#> [1] 1
#> 
#> 

# }
```
