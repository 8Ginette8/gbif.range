# Evaluate a Range Map by Cross-Validation

Rebuild a
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
model repeatedly from subsets of the occurrence data stored in a
`getRange` object and evaluate each rebuild against held-out
observations.

## Usage

``` r
cv_range(
  range_object = NULL,
  cv = "random-cv",
  nfolds = 5,
  nblocks = 2,
  backpoints = 10000,
  verbose = TRUE
)
```

## Arguments

- range_object:

  `getRange` object. Typically returned by
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md).

- cv:

  Character. String specifying the cross-validation strategy:
  `"random-cv"` or `"block-cv"`.

- nfolds:

  Numeric. Number of folds.

- nblocks:

  Numeric. Multiplier used when `cv = "block-cv"` to define the total
  number of spatial blocks as `nfolds * nblocks`.

- backpoints:

  Numeric. Number of regularly spaced background points used as
  pseudo-absences. Default is `10000`.

- verbose:

  Logical. Should fold progress messages be printed?

## Value

A data frame with one row per fold plus a `Mean` row, and the columns
`TP`, `FA`, `TA`, `FP`, `Precision`, `Sensitivity`, `Specificity`, and
`TSS`.

## Details

The function rebuilds the range map `nfolds` times. In each iteration,
one fold is reserved for evaluation and the remaining folds are used for
training.

Two strategies are available: random cross-validation and spatial block
cross-validation. The latter reduces the influence of spatial
autocorrelation by grouping nearby observations before splitting them
across folds.

Because true absences are generally unavailable, the evaluation uses a
regular grid of background points as pseudo-absences and reports
precision, sensitivity, specificity, and TSS.

## References

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-
Arroita, G., ... & Dormann, C. F. (2017). Cross-validation strategies
for data with temporal, spatial, hierarchical, or phylogenetic
structure. Ecography, 40(8), 913-929.
[doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881)

Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., &
Thuiller, W. (2021). Novel methods to correct for observer and sampling
bias in presence-only species distribution models. Global Ecology and
Biogeography, 30(11), 2312-2325.
[doi:10.1111/geb.13383](https://doi.org/10.1111/geb.13383)

## See also

[`get_range`](https://8ginette8.github.io/gbif.range/reference/get_range.md)()
to build the range map being evaluated, and
[`make_blocks`](https://8ginette8.github.io/gbif.range/reference/make_blocks.md)()
for the underlying fold-assignment logic.

## Examples

``` r
# \donttest{
###########################################
### Example plot
###########################################

# Load available ecoregions
eco.terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir(),
    format = "sf"
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
    format = "sf"
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
pt.test <- cv_range(
    range_object = range.tiger,
    cv = "block-cv",
    nfolds = 5,
    nblocks = 2
);pt.test
#> 6 variables with 5, 5, 5, 5, ... levels: 15625 function evaluations required.
#> ...fold1
#> ...fold2
#> ...fold3
#> ...fold4
#> ...fold5
#> 
#>         TP    FA     TA   FP Precision Sensitivity Specificity       TSS
#> CV1  432.0 105.0 3190.0 79.0 0.8454012   0.8044693   0.9758336 0.7803029
#> CV2  435.0 102.0 2908.0 76.0 0.8512720   0.8100559   0.9745308 0.7845867
#> CV3  256.0 281.0 2656.0 38.0 0.8707483   0.4767225   0.9858946 0.4626171
#> CV4  536.0   2.0   85.0 73.0 0.8801314   0.9962825   0.5379747 0.5342572
#> CV5  255.0 254.0  835.0 60.0 0.8095238   0.5009823   0.9329609 0.4339432
#> Mean 382.8 148.8 1934.8 65.2 0.8514153   0.7177025   0.8814389 0.5991414

###########################################
### Package manuscript plot (Fig 2a-b)
###########################################

# Root and package
root_dir <- tempdir()
if (!dir.exists(file.path(root_dir, "fig_plots"))) {
    dir.create(file.path(root_dir, "fig_plots"))
}

# -------------------------------------
# Plant
# -------------------------------------

# Preliminary
spdf.world <- terra::vect(
  system.file("extdata", "world_countries.shp", package = "gbif.range")
)
shp.lonlat <- terra::vect(
    paste0(
        system.file(package = "gbif.range"),
        "/extdata/shp_lonlat.shp"
    )
)
obs.arcto <- get_gbif("Arctostaphylos alpinus", geo = shp.lonlat, grain = 1)
#> |--------------------------------------------|
#> | Total number (all records)    :      42972 |
#> | Kept records                  :       6340 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE by default ('geo' was set)
#> 
#> ...GBIF records of Arctostaphylos alpinus: download starting...
#> ------------- #1 (100%..)               
#> Too many requests! To download GBIF occurrence data in bulk, please use occ_download().
#> Waiting [=    ] 1/5 secWaiting [==   ] 2/5 secWaiting [===  ] 3/5 secWaiting [==== ] 4/5 secWaiting [=====] 5/5 sec
#> 
#> ...Records (XY) filtering summary:
#> -----------------------------------------------
#>                     step removed remaining
#>          Grain filtering    4968      1372
#>       Duplicated records     448       924
#>          Absence records       0       924
#>          Basis selection      87       837
#>  Establishment selection       0       837
#>               Time frame       0       837
#>        Identical records       0       837
#>         Raster centroids       0       837
#> 
#> Initial records         : 6340
#> Total removed           : 5503
#> Final records (XY)      : 837
#> -----------------------------------------------
#> Final records (no XY)   : 0
rst <- terra::rast(
    paste0(
        system.file(package = "gbif.range"),
        "/extdata/rst.tif"
    )
)
my.eco <- make_ecoreg(rst, 200)
#> CLARA algorithm processing... 
#> Generating polygons... 
range.arcto <- get_range(
    occ_coord = obs.arcto,
    ecoreg = my.eco,
    ecoreg_name = "EcoRegion",
    res = 0.05,
    format = "SpatRaster"
)
#> ## Start of computation for species:  Arctostaphylos alpinus  ### 
#> 0 outlier's from 807 | proportion from total points: 0%
#> ecoregion 1  of  86 :  100 
#> ecoregion 2  of  86 :  101 
#> ecoregion 3  of  86 :  102 
#> ecoregion 4  of  86 :  103 
#> ecoregion 5  of  86 :  104 
#> ecoregion 6  of  86 :  105 
#> ecoregion 7  of  86 :  106 
#> ecoregion 8  of  86 :  107 
#> ecoregion 9  of  86 :  108 
#> ecoregion 10  of  86 :  109 
#> ecoregion 11  of  86 :  110 
#> ecoregion 12  of  86 :  111 
#> ecoregion 13  of  86 :  112 
#> ecoregion 14  of  86 :  113 
#> ecoregion 15  of  86 :  114 
#> ecoregion 16  of  86 :  115 
#> ecoregion 17  of  86 :  116 
#> ecoregion 18  of  86 :  117 
#> ecoregion 19  of  86 :  118 
#> ecoregion 20  of  86 :  119 
#> ecoregion 21  of  86 :  120 
#> ecoregion 22  of  86 :  121 
#> ecoregion 23  of  86 :  122 
#> ecoregion 24  of  86 :  123 
#> ecoregion 25  of  86 :  124 
#> ecoregion 26  of  86 :  126 
#> ecoregion 27  of  86 :  127 
#> ecoregion 28  of  86 :  128 
#> ecoregion 29  of  86 :  129 
#> ecoregion 30  of  86 :  130 
#> ecoregion 31  of  86 :  131 
#> ecoregion 32  of  86 :  133 
#> ecoregion 33  of  86 :  134 
#> ecoregion 34  of  86 :  136 
#> ecoregion 35  of  86 :  137 
#> ecoregion 36  of  86 :  138 
#> ecoregion 37  of  86 :  141 
#> ecoregion 38  of  86 :  142 
#> ecoregion 39  of  86 :  143 
#> ecoregion 40  of  86 :  144 
#> ecoregion 41  of  86 :  147 
#> ecoregion 42  of  86 :  148 
#> ecoregion 43  of  86 :  151 
#> ecoregion 44  of  86 :  189 
#> ecoregion 45  of  86 :  21 
#> ecoregion 46  of  86 :  22 
#> ecoregion 47  of  86 :  23 
#> ecoregion 48  of  86 :  25 
#> ecoregion 49  of  86 :  28 
#> ecoregion 50  of  86 :  32 
#> ecoregion 51  of  86 :  33 
#> ecoregion 52  of  86 :  36 
#> ecoregion 53  of  86 :  37 
#> ecoregion 54  of  86 :  39 
#> ecoregion 55  of  86 :  49 
#> ecoregion 56  of  86 :  52 
#> ecoregion 57  of  86 :  54 
#> ecoregion 58  of  86 :  55 
#> ecoregion 59  of  86 :  57 
#> ecoregion 60  of  86 :  58 
#> ecoregion 61  of  86 :  65 
#> ecoregion 62  of  86 :  66 
#> ecoregion 63  of  86 :  7 
#> ecoregion 64  of  86 :  70 
#> ecoregion 65  of  86 :  71 
#> ecoregion 66  of  86 :  72 
#> ecoregion 67  of  86 :  73 
#> ecoregion 68  of  86 :  74 
#> ecoregion 69  of  86 :  76 
#> ecoregion 70  of  86 :  79 
#> ecoregion 71  of  86 :  83 
#> ecoregion 72  of  86 :  84 
#> ecoregion 73  of  86 :  85 
#> ecoregion 74  of  86 :  87 
#> ecoregion 75  of  86 :  88 
#> ecoregion 76  of  86 :  89 
#> ecoregion 77  of  86 :  90 
#> ecoregion 78  of  86 :  91 
#> ecoregion 79  of  86 :  92 
#> ecoregion 80  of  86 :  93 
#> ecoregion 81  of  86 :  94 
#> ecoregion 82  of  86 :  95 
#> ecoregion 83  of  86 :  96 
#> ecoregion 84  of  86 :  97 
#> ecoregion 85  of  86 :  98 
#> ecoregion 86  of  86 :  99 
#> ## End of computation for species:  Arctostaphylos alpinus  ### 
ext.temp <- terra::ext(range.arcto$rangeOutput)
ext.temp <- c(ext.temp[1]-0.6, ext.temp[2]+0.05,
                ext.temp[3]-0.05, ext.temp[4]+0.05)

# Create pseudo-absences
    # First remove observations considered outliers in get_range()
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
#> 6 variables with 5, 5, 5, 5, ... levels: 15625 function evaluations required.
all.xy$bcv <- cv.strat
all.xy[all.xy$bcv%in%1, "col"] <- "#e41a1c"
all.xy[all.xy$bcv%in%2, "col"] <- "#377eb8"
all.xy[all.xy$bcv%in%3, "col"] <- "#4daf4a"
all.xy[all.xy$bcv%in%4, "col"] <- "#984ea3"
all.xy[all.xy$bcv%in%5, "col"] <- "#ff7f00"
all.xy[all.xy$Pres%in%1, "col"] <- substr(
  grDevices::adjustcolor(all.xy[all.xy$Pres%in%1, "col"],
    red.f = 0.5, green.f = 0.5, blue.f = 0.5), 1, 7
)

# Plot
    # Evaluate
ar.test <- cv_range(
    range_object = range.arcto,
    cv = "block-cv",
    nfolds = 5,
    nblocks = 2
)
#> 6 variables with 5, 5, 5, 5, ... levels: 15625 function evaluations required.
#> ...fold1
#> ...fold2
#> ...fold3
#> ...fold4
#> ...fold5
#> 

# First extract the world at the arcto extent and divide
    # presences/absences/outliers
world.local <- terra::crop(spdf.world, ext.temp)
world.local.ar <- terra::aggregate(world.local)
pres <- all.xy[all.xy$Pres%in%1, ]
abs <- all.xy[all.xy$Pres%in%0, ]
pres_coords <- as.data.frame(pres[, c("decimalLongitude", "decimalLatitude")])
id.in <- terra::extract(range.arcto$rangeOutput, pres_coords)
pres <- pres[!is.na(id.in[,2]), ]

    # Continue plotting
grDevices::png(
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
oldpar <- graphics::par(
    mfrow = c(1,1),
    mar = c(5,5,5,20),
    lwd = 1,
    cex = 0.5
)
terra::plot(world.local.ar, col = "#bcd1bc", axes = FALSE, lwd = 2)
graphics::points(abs, col = paste0(abs$col, "50"), pch = 16, cex = 0.7)
terra::plot(
    terra::as.polygons(range.arcto$rangeOutput),
    border = "black",
    lwd = 6,
    col = "#63636370",
    add = TRUE
)
graphics::points(pres, col = paste0(pres$col,"90"), pch = 16, cex = 1.5)
graphics::text(
    6.6,43.2,
    paste("Mean TSS =", round(tail(ar.test[,"TSS"],1),2)),
    cex = 1.5,
    font = 2
)
graphics::text(
    7.2,42.8,
    paste("Mean Precision =", round(tail(ar.test[,"Precision"],1),2)),
    cex = 1.5,
    font = 2
)
grDevices::dev.off()
#> agg_record_1ded4df2731d 
#>                       2 
graphics::par(oldpar)

# -------------------------------------
# Tiger
# -------------------------------------

# Convert the range polygon to raster
tiger_rast <- terra::rasterize(
    range.tiger$rangeOutput,
    terra::rast(
        terra::ext(range.tiger$rangeOutput),
        resolution = 0.1
    )
)

# Preliminary
ext.temp <- terra::ext(tiger_rast)
ext.temp <- c(ext.temp[1]-0.2, ext.temp[2]+0.2,
                ext.temp[3]-0.2, ext.temp[4]+0.2)

# Create pseudo-absences
    # First remove observations considered outliers in get_range()
xy.df <- range.tiger$init.args$occ_coord
r.ext <- terra::ext(tiger_rast)
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
#> 6 variables with 5, 5, 5, 5, ... levels: 15625 function evaluations required.
all.xy$bcv <- cv.strat
all.xy[all.xy$bcv%in%1,"col"] <- "#e41a1c"
all.xy[all.xy$bcv%in%2,"col"] <- "#377eb8"
all.xy[all.xy$bcv%in%3,"col"] <- "#4daf4a"
all.xy[all.xy$bcv%in%4,"col"] <- "#984ea3"
all.xy[all.xy$bcv%in%5,"col"] <- "#ff7f00"
all.xy[all.xy$Pres%in%1, "col"] <- substr(
  grDevices::adjustcolor(all.xy[all.xy$Pres%in%1, "col"],
    red.f = 0.5, green.f = 0.5, blue.f = 0.5), 1, 7
)

# Plot
    # First extract the world at the tiger extent and divide
    # presences/absences/outliers
world.local <- terra::crop(spdf.world, ext.temp)
world.local.ti <- terra::aggregate(world.local)
pres <- all.xy[all.xy$Pres%in%1, ]
abs <- all.xy[all.xy$Pres%in%0, ]
pres_coords <- as.data.frame(pres[, c("decimalLongitude", "decimalLatitude")])
id.in <- terra::extract(tiger_rast, pres_coords)
pres <- pres[!is.na(id.in[,2]),]

    # Continue plotting
grDevices::png(
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
oldpar <- graphics::par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 1, cex = 0.5)
terra::plot(world.local.ti, col = "#bcd1bc", axes = FALSE, lwd = 2)
graphics::points(abs, col = paste0(abs$col,"50"), pch = 16, cex = 0.6)
terra::plot(
    terra::as.polygons(tiger_rast),
    border = "black",
    lwd = 6,
    col = "#63636370",
    add = TRUE
)
graphics::points(pres, col = paste0(pres$col,"80"), pch = 16, cex = 1.6)
graphics::text(86,49,
    paste("Mean TSS =", round(tail(pt.test[,"TSS"],1),2)),
    cex = 1.5,
    font = 2
)
graphics::text(90.4,45.4,
    paste("Mean Precision =", round(tail(pt.test[,"Precision"],1),2)),
    cex = 1.5, font = 2
)
grDevices::dev.off()
#> agg_record_1ded4df2731d 
#>                       2 
graphics::par(oldpar)

# }
```
