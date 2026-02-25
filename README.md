# gbif.range R package

[<img align="right" width="250" height="290" src="https://github.com/8Ginette8/gbif.range/blob/main/inst/logo/logo_gbif.range.png">](https://www.gbif.org)

[![Auto-Version](https://github.com/8Ginette8/gbif.range/actions/workflows/R-Package-Auto-Version.yml/badge.svg?branch=main)](https://github.com/8Ginette8/gbif.range/actions/workflows/R-Package-Auto-Version.yml)
[![R-CMD-check](https://github.com/8Ginette8/gbif.range/actions/workflows/R-CMD-check-month-test.yml/badge.svg?branch=main)](https://github.com/8Ginette8/gbif.range/actions/workflows/R-CMD-check-month-test.yml)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Although species ranges may be obtained using expert maps (e.g., <a href="https://www.iucnredlist.org/resources/spatial-data-download">IUCN</a> and <a href="https://www.euforgen.org/species/">EUFORGEN</a>) or modeling methods, expert data remains limited in the number of available species while applying models usually need more technical expertise, as well as many species observations.

When unavailable, such information may be extracted from the Global Biodiversity Information facility (GBIF), the largest public data repository inventorying georeferenced species observations worldwide (https://www.gbif.org/). However, retrieving GBIF records at large scale in R may be tedious, if users are unaware of the limitations of the *rgbif* library.

Here we present **gbif.range**, a R library that contains automated methods to generate species range maps from scratch using in-house ecoregions shapefiles and an easy-to-use GBIF download wrapper. Finally, this library also offers a set of additional very useful tools for large GBIF datasets (generate doi, extract GBIF taxonomy, records filtering...).

_(source: globe image from the Noun Project adapted by LenaCassie-Studio)_

## Main functions

  - `get_gbif()`: improves the accessibility of the `rgbif` R package (<a href="https://cran.r-project.org/web/packages/rgbif/index.html">CRAN</a>) in
  retrieving GBIF observations of a given species (accepted and synonym names). It uses a dynamic moving windows if the given geographic extent
  contains > 100,000 observations and implements 13 post-processing options to flag and clean erroneous records based on custom functions and the
  `CoordinateCleaner` R package (<a href="https://cran.r-project.org/web/packages/CoordinateCleaner/index.html">CRAN</a>).

  - `get_range()`: estimates species ranges based on occurrence data (a `getGBIF` output or a set of coordinates) and
  <a href="https://en.wikipedia.org/wiki/Ecoregion">ecoregion</a> polygons.

  - `read_ecoreg()`: download and read available ecoregion files from different available URL sources. See also associated calls `ecoreg_list`, `get_ecoreg()` and `check_and_get_ecoreg()`.
    
  - `get_status()`: generates, based on a given species name, its IUCN red list status and a list of all scientific names
  (accepted, synonyms) found in the GBIF backbone taxonomy. Children and related doubtful names not used to download the data may also be extracted.

  - `obs_filter()`: `obs_filter()` accepts as input a `getGBIF` output (one or several species) and filter the observations according
  to a specific given grid resolution. It can retain one observation per grid pixel and/or remove observations from grid pixels that contain fewer than a specified number of records.

  - `make_tiles()`: may be used to generate a set of `SpatialExtent` and geometry arguments `POLYGON()` based on a given
  geographic extent. This function is meant to help users who want to use the `rgbif` R package and its parameter
  `geometry` that uses a `POLYGON()` argument.

  - `get_doi()`: a small wrapper of `derived_dataset()` in `rgbif` that simplifies the obtention of a general DOI
  for a set of several gbif species datasets.

  - `make_ecoreg()`: a function to create custom ecoregions based on environmental layers.

  - `evaluate_range()`: evaluation function to validate the species ranges with distribution information provided by the user.

  - `cv_range()`: cross-validation function to evaluate a `getRange` output based on its occurrence data.

## Installation

You can install the development version from GitHub with (make sure the R package `remotes` is up to date):

``` r
remotes::install_github("8Ginette8/gbif.range")
library(gbif.range)
```

## Example

### Terrestrial species

Let's download worldwide the records of *Panthera tigris* only based on true observations and literature (default):

``` r
# Download
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Plot species records
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(countries, col = "#bcbddc")
points(obs.pt[, c("decimalLongitude","decimalLatitude")], pch = 20, col = "#99340470", cex = 1.5)
```

![image](https://github.com/user-attachments/assets/2eef058d-3a0f-4b63-b053-7b1e4a94f3a2)


Note that the function did not manage to get rid of observations of most likely non-informed captive individuals (e.g., in Europe, U.S. and South Africa); see the `CoordinateCleaner` R package (<a href="https://cran.r-project.org/web/packages/CoordinateCleaner/index.html">CRAN</a>) for improved filtering. We can also retrieve the tiger **IUCN red list status**, and its scientific names (accepted and synonyms) that were used in the download with the **GBIF backbone taxonomy**. If all = TRUE, additional children and related doubtful names may also be extracted (not used in `get_gbif()`):

``` r
get_status("Panthera tigris", all = FALSE)
```

Let's now extract the terrestrial ecoregions of the world (Nature Conservancy) and generate the distributional range map of *Panthera tigris* :

``` r
# Download ecoregion and read
eco.terra <- read_ecoreg(ecoreg_name = "eco_terra", save_dir = NULL)

# Range
range.tiger <- get_range(occ_coord = obs.pt,
                        ecoreg = eco.terra,
                        ecoreg_name = "ECO_NAME",
                        degrees_outlier = 5,
                        clust_pts_outlier = 4)
```

Let's plot the result now:

``` r
terra::plot(countries, col = "#bcbddc")
terra::plot(range.tiger$rangeOutput, col = "#238b45", add = TRUE, axes = FALSE, legend = FALSE)
```

![image](https://github.com/user-attachments/assets/97f6d768-4858-42ca-89f4-9d66e2394171)


Here, default parameters were employed, however, `clust_pts_outlier` (in degrees, ~440 km here) could have been increased to remove larger isolated clusters of observations, and `degrees_outlier` (~550 km here) to account for more apart observations in the range process. Here, default parameters still allowed to remove obvious tiger observation anomalies in Europe, U.S. and South Africa.

### Available ecoregions

Although whatever shapefile may be set in `get_range()` as input, note that ecoregion shapefiles may be downloaded using the package: *eco.earth* (for terrestrial species; The Nature conservancy 2009 adapted from Olson & al. 2001), *eco.marine* (for marine species, two versions; The Nature Conservancy 2012 adapted from Spalding & al. 2007, 2012) and *eco.fresh* (for freshwater species; Abell & al. 2008). Each are available under different precision levels:
- *eco_terra* has three different levels: 'ECO_NAME', 'WWF_MHTNAM' and 'WWF_REALM2'.
- *eco_fresh* has only one: 'ECOREGION'.
- *eco_marine* and *eco_hd_marine* (very coastal-precise version) contains three distinct levels: 'ECOREGION', 'PROVINCE' and 'REALM'.

Available ecoregion files that can be downloaded with the package:
``` r
# List
ecoreg_list
```

### Custom ecoregions

Additionally, if the in-house ecoregions are too coarse for a given geographic region (e.g., for local studies) or an ecoshapefile of finer environmental details is needed, `make_ecoreg()` can be used based on spatially-informed data (e.g. climate, biodiversity) of desired resolution and extent defining the study area.

Example of 10 ecoregions based on CHELSA bioclimatic layers at 5 × 5 km resolution (Karger et al. 2017), i.e., mean annual air temperature (bio1) and annual precipitation amount (bio12) 1981–2010:

``` r
bio <- terra::rast(paste0(system.file(package = "gbif.range"), "/extdata/rst.tif"))
eco.eg <- make_ecoreg(env = bio, nclass = 10)
terra::plot(eco.eg, col = rainbow(10))
```
<img width="450" height="auto" alt="image" src="https://github.com/user-attachments/assets/3e26fe66-fbdd-4854-a24e-125a02d3928e" />

Let's further demonstrate how a custom map of ecoregions can be employed in combination with the package's main functions:

``` r
# Let's download the observations of Arctostaphylos alpinus in the European Alps:
shp.lonlat <- terra::vect(paste0(system.file(package = "gbif.range"), "/extdata/shp_lonlat.shp"))
obs.arcto <- get_gbif(sp_name = "Arctostaphylos alpinus",
                      geo = shp.lonlat,
                      grain = 1)

# Create an ecoregion layer of 200 classes, based on two environmental spatial layers:
rst <- terra::rast(paste0(system.file(package = "gbif.range"), "/extdata/rst.tif"))
my.eco <- make_ecoreg(env = rst,
                        nclass = 200)

# Create the range map based on our custom ecoregion
# (always set 'EcoRegion' as a name when using a make_ecoreg() output):
range.arcto <- get_range(occ_coord = obs.arcto,
                        ecoreg = my.eco,
                        ecoreg_name = "EcoRegion",
                        degrees_outlier = 5,
                        clust_pts_outlier = 4,
                        res = 0.05)
```

Unlike at larger-scales, we have here decreased the `get_gbif()` `grain` parameter from 100km to 1km, as keeping observations with a precision of 100km would have been too coarse to infer the approximate range distribution of the species relative to the study extent. `degrees_outlier` and `clust_pts_outlier` were here also kept defaults (~550 and 440 km, respectively), so relative to the study extent, almost no clustered or too distance observations were considered outliers.

It is also important to note that the resolution parameter (`res`) can be changed to adjust how fine the spatial output should be. This highest possible resolution will only depend on the precision of the `ecoreg` object (e.g., a range output can reach the same resolution of the rasters used to create a `make_ecoreg` object).

``` r
# Plot
alps.shp <- terra::crop(countries,terra::ext(rst))
r.arcto <- terra::mask(range.arcto$rangeOutput,alps.shp)
terra::plot(alps.shp, col = "#bcbddc")
terra::plot(r.arcto, add = TRUE, col = "darkgreen", axes = FALSE, legend = FALSE)
points(obs.arcto[, c("decimalLongitude","decimalLatitude")], pch = 20, col = "#99340470", cex = 1)
```

<img width="559" height="394" alt="image" src="https://github.com/user-attachments/assets/39928ecc-d5b8-4063-8980-cf5efae2833f" />


### Marine species

Let's reapply the same process as for Panthera tigris, but with the marine species *Delphinus delphis* (> 100'000 observations).

⚠️Notes that the download takes here longer unless the parameter `occ_samp` is used. Altough giving **less precise observational distribution**, `occ_samp` allows to extract a **subsample of *n* GBIF observations** per created tiles over the study area:

``` r
# Here the example is a sample of 1000 observations per geographic tile
obs.dd <- get_gbif("Delphinus delphis", occ_samp = 1000)

# Here the list is longer because 'all=TRUE' includes every names (even doubtful)
get_status("Delphinus delphis", all = TRUE)
```

Let's now generate three range maps of *Delphinus delphis* using the *eco.marine* as ecoregion shapefile:

``` r
# Download ecoregion and read
eco.marine <- read_ecoreg(ecoreg_name = "eco_marine", save_dir = NULL)

# Range from different levels
range.dd1 <- get_range(obs.dd, eco.marine, "ECOREGION")
range.dd2 <- get_range(obs.dd, eco.marine, "PROVINCE")
range.dd3 <- get_range(obs.dd, eco.marine, "REALM")
```

The three results are pretty similar because most of the observations are near the coast. But let's plot the first more fine result:

``` r
terra::plot(countries, col = "#bcbddc")
terra::plot(range.dd3$rangeOutput, col = "#238b45", add = TRUE, axes = FALSE, legend = FALSE)
points(obs.dd[, c("decimalLongitude","decimalLatitude")], pch = 20, col = "#99340470", cex = 1)
```

![image](https://github.com/user-attachments/assets/52d63434-f64b-4076-b8bc-d3c03c899137)

Although our result map follows the sampling pattern found in <a href="https://www.gbif.org/species/8324617">GBIF</a>, the dolphin range map might have been improved if more GBIF observations would have been extracted. Therefore, `occ_samp` must be in this case increased or removed.

## Citation

Yohann Chauvier, Oskar Hagen, Stefan Pinkert, Camille Albouy, Fabian Fopp, Philipp Brun, Patrice Descombes, Florian Altermatt, Loic Pellissier, Katalin Csilléry. gbif.range: An R package to generate ecologically-informed species range maps from occurrence data with seamless GBIF integration. Authorea. June 30, 2025.
doi: <a href="https://doi.org/10.22541/au.175130858.83083354/v1">10.22541/au.175130858.83083354/v1</a>

## References 

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity information facility API. doi: <a href="https://doi.org/10.5281/zenodo.6023735">10.5281/zenodo.6023735</a>

Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, H.P., , Kessler, M. (2017). Climatologies at high resolution for the earth’s land surface areas. Sci Data 4, 170122. doi: <a href="https://doi.org/10.1038/sdata.2017.122">10.1038/sdata.2017.122</a>

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5), 744-751. doi: <a href="https://doi.org/10.1111/2041-210X.13152">10.1111/2041-210X.13152</a>

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Link to package: <a href="https://cran.r-project.org/web/packages/terra/index.html">terra - CRAN</a>

Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. BioScience 51(11):933-938.  doi: <a href="https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2">10.1641/0006-3568(2001)051</a>

The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types, Biogeographical Realms and The Nature Conservancy Terrestrial Assessment Units. GIS layers developed by The Nature Conservancy with multiple partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986. Cambridge (UK): The Nature Conservancy. Data URL: https://geospatial.tnc.org/datasets/b1636d640ede4d6ca8f5e369f2dc368b/about

Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A. Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana, Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer Molnar, Cheri A. Recchia, James Robertson, Marine Ecoregions of the World: A Bioregionalization of Coastal and Shelf Areas, BioScience, Volume 57, Issue 7, July 2007, Pages 573–583. doi: <a href="https://doi.org/10.1641/B570707">10.1641/B570707</a>

Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012). Pelagic provinces of the world: a biogeographic classification of the world’s surface pelagic waters. Ocean & Coastal Management, 60, 19-30. doi: <a href="https://doi.org/10.1016/j.ocecoaman.2011.12.016">10.1016/j.ocecoaman.2011.12.016</a>

The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces of the World. GIS layers developed by The Nature Conservancy with multiple partners, combined from Spalding et al. (2007) and Spalding et al. (2012). Cambridge (UK): The Nature Conservancy. Data URL: http://data.unep-wcmc.org/datasets/38

Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice Kottelat, Nina Bogutskaya, Brian Coad, Nick Mandrak, Salvador Contreras Balderas, William Bussing, Melanie L. J. Stiassny, Paul Skelton, Gerald R. Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf, James Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel, Eric Wikramanayake, David Olson, Hugo L. López, Roberto E. Reis, John G. Lundberg, Mark H. Sabaj Pérez, Paulo Petry, Freshwater Ecoregions of the World: A New Map of Biogeographic Units for Freshwater Biodiversity Conservation, BioScience, Volume 58, Issue 5, May 2008, Pages 403–414. doi: <a href="https://doi.org/https://doi.org/10.1641/B580507">10.1641/B580507</a>
