# gbif.range R package

[<img align="right" width="250" height="290" src="https://github.com/8Ginette8/gbif.range/blob/main/inst/logo/logo_gbif.range.png">](https://www.gbif.org)

Status of the automatic CI R-CMD-check test

[![R-CMD-check-test](https://github.com/8Ginette8/gbif.range/actions/workflows/R-CMD-check-month-test.yml/badge.svg?branch=main&event=push)](https://github.com/8Ginette8/gbif.range/actions/workflows/R-CMD-check-month-test.yml)


Although species ranges may be obtained using expert maps (e.g., <a href="https://www.iucnredlist.org/resources/spatial-data-download">IUCN</a> and <a href="https://www.euforgen.org/species/">EUFORGEN</a>) or modeling methods, expert data remains limited in the number of available species while applying models usually need more technical expertise, as well as many species observations.

When unavailable, such information may be extracted from the Global Biodiversity Information facility (GBIF), the largest public data repository inventorying georeferenced species observations worldwide (https://www.gbif.org/). However, retrieving GBIF records at large scale in R may be tedious, if users are unaware of the limitations of the *rgbif* library.

Here we present **gbif.range**, a R library that contains automated methods to generate species range maps from scratch using in-house ecoregions shapefiles and an easy-to-use GBIF download wrapper. Finally, this library also offers a set of additional very useful tools for large GBIF datasets (generate doi, extract GBIF taxonomy, records filtering...).

_(source: globe image from the Noun Project adapted by LenaCassie-Studio)_

## Main functions

  - *get_gbif()*: improves the accessibility of the *rgbif* R package (<a href="https://cran.r-project.org/web/packages/rgbif/index.html">CRAN</a>) in
  retrieving GBIF observations of a given species (accepted and synonym names). It uses a dynamic moving   windows if the given geographic extent
  contains > 100,000 observations and implements 13 post-processing options to flag and clean erroneous records based on custom functions and the
  *CoordinateCleaner* R package (<a href="https://cran.r-project.org/web/packages/CoordinateCleaner/index.html">CRAN</a>).

  - *get_range()*: estimates species ranges based on occurrence data (a *get_gbif* output or a set of coordinates) and
  <a href="https://en.wikipedia.org/wiki/Ecoregion">ecoregion</a> polygons.

  - *read_bioreg()*: download and read available ecoregion files from different available URL sources. See also associated calls *bioreg_list*, *get_bioreg()* and *check_and_get_bioreg()*.
    
  - *get_status()*: generates, based on a given species name, its IUCN red list status and a list of all scientific names
  (accepted, synonyms) found in the GBIF backbone taxonomy. Children and related doubtful names not used to download the data may also be extracted.

  - *obs_filter()*: *obs_filter()* accepts as input a *get_gbif()* output (one or several species) and filter the observations according
  to a specific given grid resolution (one observation per pixel grid kept). This function allows the user to refine the density of GBIF
  observations according to a defined analysis/study's resolution.

  - *make_tiles()*: may be used to generate a set of *SpatialExtent* and geometry arguments POLYGON() based on a given
  geographic extent. This function is meant to help users who want to use the *rgbif* R package and its parameter
  *geometry* that uses a POLYGON() argument.

  - *get_doi()*: a small wrapper of *derived_dataset()* in *rgbif* that simplifies the obtention of a general DOI
  for a set of several gbif species datasets.

  - *make_ecoregion()*: a function to create custom ecoregions based on environmental layers.

  - *evaluate_range()*: evaluation function to validate the species ranges with distribution information provided by the user.

  - *cv_range()*: cross-validation function to evaluate a *get_range()* output based on its occurrence data.

## Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("8Ginette8/gbif.range")
library(gbif.range)
library(terra)
library(rnaturalearth)
```

## Example

### Terrestrial species

Let's download worldwide the records of *Panthera tigris* only based on true observations and litterature information:

``` r
# Download
obs_pt = get_gbif(sp_name = "Panthera tigris")

# Plot species records
countries = ne_countries(type = "countries",returnclass = "sv")
plot(countries,col = "#bcbddc")
points(obs_pt[,c("decimalLongitude","decimalLatitude")],pch=20,col="#99340470",cex=1.5)
```

![image](https://github.com/user-attachments/assets/2af40d4a-f1e5-47ba-a9ad-9e6da32a9df4)


Note that the function did not manage to get rid of observations found in the US and germany (observations from zoos most likely). We can also retrieve the tiger **IUCN red list status**, and its scientific names (accepted and synonyms) that were used in the download with the **GBIF backbone taxonomy**. If all = TRUE, additonal children and related doubtful names may also be extracted (not used in *get_gbif()*):

``` r
get_status("Panthera tigris",all=FALSE)
```

Let's now extract the terrestrial ecoregions of the world (Nature Conservancy) and generate the distributional range map of *Panthera tigris* :

``` r
# Download ecoregion and read
eco_terra = read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# Range
range_tiger = get_range(occ_coord = obs_pt,
                        bioreg = eco_terra,
                        bioreg_name = "ECO_NAME",
                        degrees_outlier = 6,
                        clustered_points_outlier = 4)
```

Let's plot the result now:

``` r
plot(countries,col = "#bcbddc")
plot(range_tiger$range_output,col = "#238b45",add = TRUE,axes = FALSE,legend = FALSE)
```

![image](https://github.com/user-attachments/assets/6cf6c8c3-7e0f-4754-9806-05821091173c)

Here, *clustered_points_outlier = 4* was employed to remove US and European observations of *Panthera tigris* from the range process, and *degrees_outlier* slightly increased to *6* to account for more appart observations in the range process. Note that buffer and filtering parameters can be be set in *get_range* and should be carefully explored before any definite range map generation.

### Available ecoregions

Although whatever shapefile may be set in *get_range()* as input, note that ecoregion shapefiles may be dowload using the package: *eco.earh* (for terrestrial species; The Nature conservancy 2009 adapted from Olson & al. 2001), *eco.marine* (for marine species, two versions; The Nature Conservancy 2012 adapted from Spalding & al. 2007, 2012) and *eco.fresh* (for freshwater species; Abell & al. 2008). Each are available under different precision levels:
- *eco_terra* has three different levels: 'ECO_NAME', 'WWF_MHTNAM' and 'WWF_REALM2'.
- *eco_fresh* has only one: 'ECOREGION'.
- *eco_marine* and *eco_hd_marine* (very coastal-precise version) contains three distinct levels: 'ECOREGION', 'PROVINCE' and 'REALM'.

Available ecoregion files that can be downloaded with the package:
``` r
# List
bioreg_list
```

### Custom ecoregions

Additonally, if the in-house ecoregions are too coarse for a given geographic region (e.g., for local studies) or an ecoshapefile of finer environmental details is needed, *make_ecoregion()* can be used based on spatially-informed environment (e.g. climate) of desired resolution and extent defining the study area; example:

``` r
# Let's download the observations of Arctostaphylos alpinus in the European Alps:
shp_lonlat = vect(paste0(system.file(package = "gbif.range"),"/extdata/shp_lonlat.shp"))
obs_arcto = get_gbif("Arctostaphylos alpinus",geo=shp_lonlat)

# Create an ecoregion layer of 200 classes, based on two environmental spatial layers:
rst = rast(paste0(system.file(package = "gbif.range"),"/extdata/rst.tif"))
my_eco = make_ecoregion(rst,200)

# Create the range map based on our custom ecoregion
# (always set 'EcoRegion' as a name when using a make_ecoregion() output):
range_arcto = get_range(occ_coord = obs_arcto,
                        bioreg = my_eco,
                        bioreg_name = "EcoRegion",
                        res = 20,
                        degrees_outlier = 5,
                        clustered_points_outlier = 3,
                        buffer_width_point = 4, 
                        buffer_increment_point_line = 0.5, 
                        buffer_width_polygon = 0.1)
```

Here we adapted the extra-parameters to the extent of the study area, e.g., (i) consider points as outliers (a maximum group of three points) if this bunch is away > 555km (1° ~ 111km) from the other cluster points and (ii) apply a buffer of ~10km around the drawn polygons. ⚠️It is also important to note that the resolution parameter ('res') can be changed to adjust how fine the spatial output should be. This highest possible resolution will only depend on the precision of the *bioreg* object (e.g., a range output can reach the same resolution of the rasters used to create a *make_ecoregion* object).

``` r
# Plot
plot(crop(countries,ext(rst)),col = "#bcbddc")
plot(range_arcto$range_output,add = TRUE,col = "darkgreen",axes = FALSE,legend = FALSE)
points(obs_arcto[,c("decimalLongitude","decimalLatitude")],pch = 20,col = "#99340470",cex=1)
```

![image](https://github.com/user-attachments/assets/832e6d57-f7cb-402a-985e-fdf05d4f96aa)

### Marine species

Let's reapply the same process as for Panthera tigris, but with the marine species *Delphinus delphis* (> 100'000 observations).

⚠️Notes that the download takes here longer unless the parameter *occ_samp* is used. Altough giving **less precise observational distribution**, *occ_samp* allows to extract a **subsample of *n* GBIF observations** per created tiles over the study area:

``` r
obs_dd = get_gbif("Delphinus delphis",occ_samp = 1000) # Here the example is a sample of 1000 observations per geographic tile
get_status("Delphinus delphis",all = TRUE) # Here the list is longer because 'all=TRUE' includes every names (even doubtful)
```

Let's now generate three range maps of *Delphinus delphis* using the *eco.marine* as ecoregion shapefile:

``` r
# Download ecoregion and read
eco_marine = read_bioreg(bioreg_name = "eco_marine", save_dir = NULL)

# Range from different levels
range_dd1 = get_range(obs_dd,eco_marine,"ECOREGION")
range_dd2 = get_range(obs_dd,eco_marine,"PROVINCE")
range_dd3 = get_range(obs_dd,eco_marine,"REALM")
```

The three results are pretty similar because most of the observations are near the coast. But let's plot the first more fine result:

``` r
plot(countries,col="#bcbddc")
plot(range_dd3$range_output,col = "#238b45",add = TRUE,axes = FALSE,legend = FALSE)
points(obs_dd[,c("decimalLongitude","decimalLatitude")],pch = 20,col = "#99340470",cex = 1)
```

<img width=80% height=80% src="https://github.com/8Ginette8/gbif.range/assets/43674773/a84c5dcf-f2c7-4722-b2ed-d13502d45eb1">

Althought our result map follows the sampling pattern found in <a href="https://www.gbif.org/species/8324617">GBIF</a>, the dolphin range map might have been improved if more GBIF observations woud have been extracted. Therefore, *occ_samp* must be in this case increased or removed.

⚠️Finally, also note that in case of too many records, *get_range* can be used with a **subsample of species observations** to ensure a **faster polygon process and/or to overcome potential RAM crashes**.

## Citation
Yohann Chauvier; Oskar Hagen; Camille Albouy; Patrice Descombes; Fabian Fopp; Michael P. Nobis; Philipp Brun; Lisha Lyu; Loïc Pellissier; Katalin Csilléry (2022). gbif.range - An R package to generate species range maps based on ecoregions and a user-friendly GBIF wrapper. EnviDat. doi: <a href="https://www.envidat.ch/#/metadata/gbif-range-r">10.16904/envidat.352</a>

## References

Oskar Hagen, Lisa Vaterlaus, Camille Albouy, Andrew Brown, Flurin Leugger, Renske E. Onstein, Charles Novaes de Santana, Christopher R. Scotese, Loïc Pellissier. (2019) Mountain building, climate cooling and the richness of cold-adapted plants in the Northern Hemisphere. Journal of Biogeography. doi: <a href="https://doi.org/10.1111/jbi.13653">10.1111/jbi.13653</a>

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. doi: <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>

Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J. S., ... & Pellissier, L. (2022). An integrated high‐resolution mapping shows congruent biodiversity patterns of Fagales and Pinales. New Phytologist, 235(2). doi: <a href="https://doi.org/10.1111/nph.18158">10.1111/nph.18158</a>

Pinkert, S., Sica, Y. V., Winner, K., & Jetz, W. (2023). The potential of ecoregional range maps for boosting taxonomic coverage in ecology and conservation. Ecography, 12, e06794. doi: <a href="https://doi.org/10.1111/ecog.06794">10.1111/ecog.06794</a>

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity information facility API. doi: <a href="https://doi.org/10.5281/zenodo.6023735">10.5281/zenodo.6023735</a>

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5), 744-751. doi: <a href="https://doi.org/10.1111/2041-210X.13152">10.1111/2041-210X.13152</a>

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Link to package: <a href="https://cran.r-project.org/web/packages/terra/index.html">terra - CRAN</a>

Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.  doi: <a href="https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2">10.1641/0006-3568(2001)051</a>

The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types, Biogeographical Realms and The Nature Conservancy Terrestrial Assessment Units. GIS layers developed by The Nature Conservancy with multiple partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986. Cambridge (UK): The Nature Conservancy. Data URL: https://geospatial.tnc.org/datasets/b1636d640ede4d6ca8f5e369f2dc368b/about

Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A. Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana, Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer Molnar, Cheri A. Recchia, James Robertson, Marine Ecoregions of the World: A Bioregionalization of Coastal and Shelf Areas, BioScience, Volume 57, Issue 7, July 2007, Pages 573–583. doi: <a href="https://doi.org/10.1641/B570707">10.1641/B570707</a>

Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012). Pelagic provinces of the world: a biogeographic classification of the world’s surface pelagic waters. Ocean & Coastal Management, 60, 19-30. doi: <a href="https://doi.org/10.1016/j.ocecoaman.2011.12.016">10.1016/j.ocecoaman.2011.12.016</a>

The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces of the World. GIS layers developed by The Nature Conservancy with multiple partners, combined from Spalding et al. (2007) and Spalding et al. (2012). Cambridge (UK): The Nature Conservancy. Data URL: http://data.unep-wcmc.org/datasets/38

Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice Kottelat, Nina Bogutskaya, Brian Coad, Nick Mandrak, Salvador Contreras Balderas, William Bussing, Melanie L. J. Stiassny, Paul Skelton, Gerald R. Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf, James Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel, Eric Wikramanayake, David Olson, Hugo L. López, Roberto E. Reis, John G. Lundberg, Mark H. Sabaj Pérez, Paulo Petry, Freshwater Ecoregions of the World: A New Map of Biogeographic Units for Freshwater Biodiversity Conservation, BioScience, Volume 58, Issue 5, May 2008, Pages 403–414. doi: <a href="https://doi.org/https://doi.org/10.1641/B580507">10.1641/B580507</a>
