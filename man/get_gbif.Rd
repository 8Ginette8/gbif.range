% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gbif.R
\name{get_gbif}
\alias{get_gbif}
\title{Massively download and filter GBIF observations for sound spatial analyses}
\usage{
get_gbif(
  sp_name = NULL,
  search = TRUE,
  rank = NULL,
  phylum = NULL,
  class = NULL,
  order = NULL,
  family = NULL,
  conf_match = 80,
  geo = NULL,
  grain = 100,
  duplicates = FALSE,
  absences = FALSE,
  no_xy = FALSE,
  basis = c("OBSERVATION", "HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OCCURRENCE",
    "MATERIAL_CITATION", "MATERIAL_SAMPLE", "LITERATURE"),
  establishment = c("native", "casual", "released", "reproducing", "established",
    "colonising", "invasive", "widespreadInvasive"),
  add_infos = NULL,
  time_period = c(1000, 3000),
  identic_xy = FALSE,
  wConverted_xy = TRUE,
  centroids = FALSE,
  ntries = 10,
  error_skip = TRUE,
  occ_samp = 99000,
  ...
)
}
\arguments{
\item{sp_name}{Character. Species name from which the user wants to retrieve all existing GBIF names
with associated taxonomy and IUCN status.}

\item{search}{Logical. If TRUE (default), the function will strictly look for the most relevant result, based
on a list of names given by rgbif (only species, subspecies and variety allowed here), and give an error
if name matching was impeded by synonym duplicates. If FALSE, the function will simply pick the first most
relevant name from the list (higher taxa level than species allowed here). Also, unlike search=TRUE,
fuzzy search (~approximative name match) is here allowed, and the 'rank', phylum', 'class', order'
and 'family' parameters are optionally used only if no convincing name match is found.}

\item{rank}{Character. "SPECIES", "SUBSPECIES" or "VARIETY". If NULL (default), the order of priority
is (1) species, (2) subspecies and (3) variety unless "subsp." or "var." is found in 'sp_name'.}

\item{phylum}{Character. Optional. What is the species' Phylum? Adds a criteria to deal with alternative
name matches and select the right synonym. Available options are the GBIF Phylums
(listed per Kingdom --> https://www.gbif.org/species/1). If search = FALSE, used only if no direct match
is found.}

\item{class}{Character. Optional. What is the species' Class? Same as above but at the finer class level.
Available options are the GBIF Classes (same url). If search = FALSE, used only if no direct match
is found.}

\item{order}{Character. Optional. What is the species' Order? Same as above but at the finer order level.
Available options are the GBIF Orders (same url). If search = FALSE, used only if no direct match
is found.}

\item{family}{Character. Optional. What is the species' Family? Same as above but at the finer family level.
Available options are the GBIF Orders (same url). If search = FALSE, used only if no direct match
is found.}

\item{conf_match}{Numeric from 0 to 100. Determine the confidence threshold of match
of 'sp_name' with the GBIF backbone taxonomy. Default is 90.}

\item{geo}{Object of class Extent, SpatExtent, SpatialPolygon, SpatialPolygonDataframe,
or SpatVector (WGS84) to define the study's area extent. Default is NULL, i.e., the whole globe.}

\item{grain}{Numeric. Specify in kilometers the study resolution. Default is 100. Used to filter gbif records
according to their (1) spatial uncertainties and (2) number of coordinate decimals. Records with no information
on coordinate uncertainties (column coordinateUncertaintyInMeters') are be kept by default. See details.}

\item{duplicates}{Logical. Should duplicated records be kept? Default is FALSE.}

\item{absences}{Logical. Should absence records be kept? Default is FALSE.}

\item{no_xy}{Logical. Only records with coordinates are downloaded. Default is FALSE.
If TRUE, records with no coordinates are also downloaded.}

\item{basis}{Character. Which basis of records should be selected?
Available (old and new) are 'OBSERVATION', 'HUMAN_OBSERVATION', 'MACHINE_OBSERVATION',
'MATERIAL_CITATION', MATERIAL_SAMPLE', 'PRESERVED_SPECIMEN', 'FOSSIL_SPECIMEN',
'LIVING_SPECIMEN', 'LITERATURE', UNKNOWN' and 'OCCURRENCE'. Default setting removes
specimens and unknown observations.
Description may be found here: https://docs.gbif.org/course-data-use/en/basis-of-record.html, 
https://gbif.github.io/parsers/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html.}

\item{establishment}{Character. Is the individual native, captive or else? Defaut is
native, casual, released, reproducing, established, colonising and absence of information.
See https://dwc.tdwg.org/list/#dwc_degreeOfEstablishment for other managed establishments:
managed, captive, cultivated, released, unestablished and failing}

\item{add_infos}{Character. Infos that may be added to the default output information.
List of IDs may be found at: https://www.gbif.org/developer/occurrence.
Default IDs contain 'taxonKey', 'scientificName', 'acceptedTaxonKey',
'acceptedScientificName', 'individualCount', 'decimalLatitude', 'decimalLongitude',
'basisOfRecord', 'coordinateUncertaintyInMeters', 'countryCode', 'country', 'year', 'datasetKey', 
'institutionCode', 'publishingOrgKey', 'taxonomicStatus', 'taxonRank' and 'degreeOfEstablishment'.}

\item{time_period}{Numerical vector. Observations will be downloaded according to the chosen
year range. Default is c(1000,3000). Observations with year = NA are kept by default.}

\item{identic_xy}{Logical. Should records with identical xy be kept? Default is FALSE.}

\item{wConverted_xy}{Logical. Should incorrectly converted lon/lat be kept? Default is TRUE.
Otherwise, implements an approximate version of cd_ddmm() from the 'CoordinateCleaner'
R package. See the package for more advanced options.}

\item{centroids}{Logical. Should species records from raster centroids be kept? Default is TRUE.
Uses cd_round() from 'CoordinateCleaner' R package.}

\item{ntries}{Numeric. In case of failure from GBIF server or within the rgbif package, how many
download attempts should the function request? Default is '10' with a 2 seconds interval
between tries. If the attempts failed, an empty data.frame is return by default.}

\item{error_skip}{Logical. Should the search process continues if ntries failed ?}

\item{occ_samp}{Numeric. Determine how many GBIF occurrences will be sampled per geographic
tiles of the fragmented study area. Default is the maximum number of GBIF observations found
in a tile (i.e. ~100'000 records). A lower number may be set (<99'000) if the user only wants
a sample of the species GBIF observations, hence increasing the download process and the
generation of its range map if get_range() is employed afterwards.}

\item{...}{Additonnal parameters for the function cd_round() of CoordinateCleaner.}
}
\value{
Object of class data.frame with requested GBIF information. Although the function
works accurately, error outputs might still occur depending on the 'sp_name' used.
Therefore, default information detailed in 'add_infos' is stored so that sanity checks
may still be applied afterwards. Although crucial preliminary checks of species records
are done by the function, additional post exploration with the CoordinateCleaner R
package is still highly recommended.
}
\description{
Implement an user-friendly workflow to download and clean gbif taxa observations.
The function uses the rgbif R package but (1) implements the same search result 
found if www.gbif.org is employed i.e., based on the input taxa name, all species
records related to its accepted name and synonyms are extracted. The function
also (2) bypasses the rgbif hard limit on the number of records (100'000 max).
For this purpose, a dynamic moving window is created and used across the geographic
extent defined by the user. This window automatically fragments the specified
study area in succesive tiles of different sizes, until all tiles include < 100'000
observations. The function also (3) automatically applies a post-filtering of
observations based on the chosen resolution of the study/analysis and by partly
employing the CoordinateCleaner R package. Filtering options may be chosen and
involve several choices: study's extent, removal of duplicates, removal of absences,
basis of records selection, removal of invalid/uncertain xy coordinates (WGS84), time
period selection and removal of raster centroids. By default, the argument
hasGeospatialIssue in occ_data() (implemented rgbif function) is set to FALSE.
}
\details{
Argument `grain` used for two distinct gbif records filtering. (1) Records filtering
according to gbif 'coordinateUncertaintyInMeters'; every records uncertainty > grain/2
are removed. Note: Records with no information on coordinate uncertainties are kept by
default. (2) Records filtering according to the number of longitude/latitude decimals;
if 110km < grain <= 11km, lon/lat with >= 1 decimal are kept, if 11km < grain <= 1100m,
lon/lat with >= 2 decimals kept; if 1100m < grain <= 110m, lon/lat with >= 3 decimals
are kept; if 110m < grain <= 11m, lon/lat with >= 4 decimals are kept;
if 11m < grain <= 1.1m, lon/lat with >= 5 decimals are kept etc...
}
\examples{
# Downloading worldwide the observations of Panthera tigris
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Plot
countries <- rnaturalearth::ne_countries(type = "countries", returnclass = "sv")
terra::plot(countries, col = "#bcbddc")
points(obs.pt[,c("decimalLongitude","decimalLatitude")], pch = 20, col = "#238b4550", cex = 4)

\dontrun{
# Downloading worldwide the observations of Ailuropoda melanoleuca (with a 100km grain, after 1990
# and by keeping duplicates and by adding the name of the person who collected the panda records)
obs.am <- get_gbif(sp_name = "Ailuropoda melanoleuca",
                   grain = 100 ,
                   duplicates = TRUE,
                   time_period = c(1990,3000),
                   add_infos = c("recordedBy","issue"))

# Plot
terra::plot(countries, col = "#bcbddc")
graphics::points(obs.am[,c("decimalLongitude","decimalLatitude")],
   pch = 20, col = "#238b4550", cex = 4)

# Downloading worlwide the observations of Phascolarctos cinereus (with a 1km grain, after 1990,
# and keeping raster centroids)
obs.pc <- get_gbif(sp_name = "Phascolarctos cinereus",
                   grain = 1,
                   time_period = c(1990,3000),
                   centroids = TRUE)
}
}
\references{
Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann,
N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the
European Alps. Ecological monographs, 91(2), e01433. 10.1002/ecm.1433

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
information facility API. 10.5281/zenodo.6023735

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli,
A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection
databases. Methods in Ecology and Evolution, 10(5), 744-751. 10.1111/2041-210X.13152

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Terra - CRAN
}
\seealso{
The (1) rgbif and (2) CoordinateCelaner packages for additional and more general
approaches on (1) downloading GBIF observations and (2) post-filtering those.
}
