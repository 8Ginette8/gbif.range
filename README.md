# gbif.range

[<img align="right" width="180" height="180" src="https://user-images.githubusercontent.com/43674773/198343094-822ff0f1-72de-4b0d-b293-ced27de96190.png">](https://www.gbif.org)

Although species range may be obtained using expert maps (IUCN - https://www.iucnredlist.org/resources/spatial-data-download, EUFORGEN - https://www.euforgen.org/species/) or modeling methods, expert data is often species-limited and statistical models need more technical expertise as well as many species observations. When unavailable, such information may be extracted from the Global Biodiversity Information facility (GBIF), the largest public data repository inventorying georeferenced species observations worldwide (https://www.gbif.org/).

However, retrieving GBIF records at large scale in R may be tedious if users are unaware of the specific set of functions and parameters that need to be employed in the *rgbif* library, and because of its existing limitation on the number of downloaded records (<100'000) if no data request is made. Here we present **gbif.range**, a R library that contains automated methods to generate species range maps from scratch using in-house ecoregions shapefiles and an easy-to-use GBIF download wrapper. Finally, this library also offers a set of additional very useful parameters and functions for large GBIF datasets (generate doi, extract GBIF taxonomy, records filtering...).

_(source: globe image adapted from Akhil Komath from the Noun Project)_

## Description
### GBIF wrapper

One the one hand, *get_gbif()* is wrapper that allows the whole observations of a given species name (accepted and synonym names) to be automatically retrieved, and improves the data accessibility of the *rgbif* R package (https://cran.r-project.org/web/packages/rgbif/index.html). The user download hard limit of *rgbif* is a maximum of 100,000 of species observations in one go if the easy-to-use interactive functions *occ_search()* and *occ_data()* are used (i.e., if no official download request is made with *occ_download()*, https://www.gbif.org/developer/occurrence). This impends the fast accessibility to the GBIF database when large observational datasets for many species and regions of the world are needed, specifically in scientific fields related to macroecology, modelling or satellite imagery. *get_gbif()* therefore bypasses this limit by intuitively using geographic parameters from *occ_data()* in *rgbif* function and the *terra* R package and adopting a dynamic moving window process allowing the user's study area of interest to be automatically fragmented in several tiles that always include < 100,000 observations.

On the other hand, *get_gbif()* also implements easy-to-use preliminary filtering options implemented during the download so that the user saves some post-processing time in data cleaning. 13 filters are available. Two already are set by default in *get_gbif()* (hasCoordinate = TRUE, hasGeospatialIssue=FALSE) and 11 can be chosen independently, including two that are based on the *CoordinateCleaner* R package (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html). It is important to note that, although a strong records filtering may be undertaken with *get_gbif()*, *CoordinateCleaner* includes a larger variety of options that should be checked and applied on *get_gbif* outputs.

### Range function

*get_range()* estimates species ranges based on occurrence data (a *get_gbif* output or a set of coordinates) and bioregions. It first deletes outliers from the observation dataset and then creates a polygon (convex hull) with a user specified buffer around all the observations of one bioregion. If there is only one observation in a bioregion, a buffer around this point will be created. If all points in a bioregion are on a line, the function will also create a buffer around these points, however, the buffer size increases with the number of points in the line. More details to come...

###  Additonal functions

*gbif.range* also includes a set of additional functions meant to be a nice supplement of the features and data that offer *get_gbif* and *get_range*:
  - *get_taxonomy*: Generates, based on a given species name, a list of all its scientific names
  (accepted, synonyms) found in the GBIF backbone taxonomy to download the data. Children and related
  doubtful names not used to download the data may also be extracted. The function allows therefore taxonomy
  correspondency to be made between different species and sub-species to potentially merge their records,
  but also permits efficient ways of linking external data of a species which is named differently across databases.
  - *obs_filter*: Whereas the 'grain' parameter in *get_gbif()* allows GBIF observations to be filtered
  according to a certain spatial precision, *obs_filter()* accepts as input a *get_gbif()* output (one or
  several species) and filter the observations according to a specific given grid resolution (one observation
  per pixel grid kept). This function allows the user to refine the density of GBIF observations according to
  a defined analysis/study's resolution.
  - *make_tiles*: May be used to generate a set of *SpatialExtent* and geometry arguments POLYGON() based on a given
  geographic extent. This function is meant to help users who want to use the *rgbif* R package and its parameter
  *geometry* that uses a POLYGON() argument.
  - *get_doi*: A small wrapper of *derived_dataset()* in *rgbif* that simplifies the obtention of a general DOI
  for a set of several gbif species datasets.

## Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("8Ginette8/gbif.range")
```

## Example

Let's download worldwide the records of Panthera tigris only based on true observations:

``` r
# Download
obs.pt = get_gbif("Panthera tigris",basis=c("OBSERVATION","HUMAN_OBSERVATION"))

# Plot species records
library(maptools)
data(wrld_simpl)
plot(wrld_simpl)
points(obs.pt[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=4)
```

Note that an additional filtering needs here to be done as one observation is found in the US. A lot of tigers are being captive in this country hence the recorded observation. Therefore *CoordinateCleaner* functions should here be considered thereafter.

We can also retrieve all the tiger scientific names (accepted and synonyms) that were used in the download with the GBIF backbone taxonomy. If all = TRUE, additonal children and related doubtful names may also be extracted (not used in *wsl_gbif()*):

``` r
get_taxonomy("Panthera tigris",all=FALSE)
```

Same may be done with Delphinus delphis (a species with > 100'00 observations)

``` r
obs.dd = get_gbif("Delphinus delphis")
get_taxonomy("Delphinus delphis",all=TRUE) # Here the list is longer because 'all=TRUE' includes every names (even doubtful)
```

Let's now generate the species range map of Panthera tigris. Although whatever shapefile may be set as input, note that three ecoregion shapefiles are already included in the library: *eco.earh* (for terrestrial species; Nature conservancy version adapted from Olson & al. 2001), *eco.marine* (for coastal and reef species; Spalding & al. 2007) and *eco.fresh* (for freshwater species; Abell & al. 2008). For deep ocean/sea species, *eco.earth* may be used, but the polygon estimates will only be geographic. Each ecoregion shapefile has one or more categories, which describe more or less precisely the ecoregion global distribution. For example, *eco.earth* has three different levels: ECO_name, WWF_MHTNAM, WWF_REALM (more to less detailed). 

``` r
names(eco.earth)
```

Which level should you pick depends on your questions and the species ecology you want to explore. Here, we choose *eco.earth* since Panthera tigris is of course a terrestrial species, and the very detailed 'ECO_NAME' as ecoregion name:

``` r
range.tiger = get_range("Panthera tigris",obs.pt,eco.earth,"ECO_NAME")
```

Other examples may be found in the R documentation *gbif.range.pdf*..

## Citation
Yohann Chauvier; Patrice Descombes; Oskar Hagen; Camille Albouy; Fabian Fopp; Michael P. Nobis; Philipp Brun; Lisha Lyu; Katalin Csilléry; Loïc Pellissier (2022). gbif.range - A R package to generate species range maps based on ecoregions and an user-friendly GBIF wrapper. EnviDat. doi: <a href="https://www.envidat.ch/#/metadata/gbif-range-r">10.16904/envidat.352.</a>

## References

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. doi: <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>

Lyu, L., Leugger, F., Hagen, O., Fopp, F., Boschman, L. M., Strijk, J. S., ... & Pellissier, L. (2022). An integrated high‐resolution mapping shows congruent biodiversity patterns of Fagales and Pinales. New Phytologist, 235(2), doi: <a href="https://doi.org/10.1111/nph.18158">10.1111/nph.18158</a>

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity information facility API. doi: <a href="https://doi.org/10.5281/zenodo.6023735">10.5281/zenodo.6023735</a>

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5), 744-751. doi: <a href="https://doi.org/10.1111/2041-210X.13152">10.1111/2041-210X.13152</a>

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Link to package: <a href="https://cran.r-project.org/web/packages/terra/index.html">Terra - CRAN</a>

Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D., Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map of life on Earth. Bioscience 51(11):933-938.  doi: <a href="https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2">10.1641/0006-3568(2001)051</a>

Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A. Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana, Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer Molnar, Cheri A. Recchia, James Robertson, Marine Ecoregions of the World: A Bioregionalization of Coastal and Shelf Areas, BioScience, Volume 57, Issue 7, July 2007, Pages 573–583. doi: <a href="https://doi.org/10.1641/B570707">10.1641/B570707</a>

Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice Kottelat, Nina Bogutskaya, Brian Coad, Nick Mandrak, Salvador Contreras Balderas, William Bussing, Melanie L. J. Stiassny, Paul Skelton, Gerald R. Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf, James Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel, Eric Wikramanayake, David Olson, Hugo L. López, Roberto E. Reis, John G. Lundberg, Mark H. Sabaj Pérez, Paulo Petry, Freshwater Ecoregions of the World: A New Map of Biogeographic Units for Freshwater Biodiversity Conservation, BioScience, Volume 58, Issue 5, May 2008, Pages 403–414. doi: <a href="https://doi.org/https://doi.org/10.1641/B580507">10.1641/B580507</a>
