# Download and Filter GBIF Occurrences for Spatial Analyses

Downloads occurrence records from GBIF for a focal taxon, resolving the
input name against the GBIF backbone taxonomy and querying by accepted
taxon key — thereby capturing records linked to synonyms and
infra-specific taxa (subspecies, varieties) in the same way as the GBIF
website. Large requests are automatically split into spatial tiles, and
a series of post-download filters can be applied to clean records for
spatial analyses and range mapping.

## Usage

``` r
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
  has_xy = TRUE,
  spatial_issue = FALSE,
  grain = 100,
  duplicates = FALSE,
  absences = FALSE,
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
  occ_samp = 10000,
  should_use_occ_download = FALSE,
  occ_download_user = NULL,
  occ_download_pwd = NULL,
  occ_download_email = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- sp_name:

  Character. String with the species name. Scientific names at
  genus-species level are expected; fuzzy matching is available when
  `search = FALSE`.

- search:

  Logical. If `TRUE` (default), use a strict GBIF backbone search and
  keep only species-, subspecies-, or variety-level matches. If `FALSE`,
  use a more permissive search and optionally rely on `rank`, `phylum`,
  `class`, `order`, and `family` to resolve ambiguous matches.

- rank:

  Character. String giving the preferred rank to keep: `"SPECIES"`,
  `"SUBSPECIES"`, or `"VARIETY"`. When `NULL`, rank priority is inferred
  from `sp_name`.

- phylum:

  Optional character. Phylum used to disambiguate alternative GBIF
  matches. Particularly useful for hemihomonyms.

- class:

  Optional character. Class used to disambiguate alternative GBIF
  matches.

- order:

  Optional character. Order used to disambiguate alternative GBIF
  matches.

- family:

  Optional character. Family used to disambiguate alternative GBIF
  matches.

- conf_match:

  Numeric. Confidence threshold between 0 and 100 for the GBIF backbone
  match. Default is `80`.

- geo:

  Spatial object. Used to restrict the query extent. Accepted classes
  are `Extent`, `SpatExtent`, `SpatialPolygon`,
  `SpatialPolygonDataFrame`, `SpatVector`, and `sf`. The default `NULL`
  queries the whole globe.

- has_xy:

  Logical. If `TRUE` (default), keep only records with coordinates. If
  `FALSE`, keep only records without coordinates. If `NULL`, keep all
  records.

- spatial_issue:

  Logical. If `FALSE` (default), keep only records without geospatial
  issues. If `TRUE`, keep only records with geospatial issues. If
  `NULL`, keep all records.

- grain:

  Numeric. Study grain in kilometers. Default is `100`. The value is
  used to filter records by `coordinateUncertaintyInMeters` and by the
  number of reported coordinate decimals.

- duplicates:

  Logical. Should duplicate records be kept? Default is `FALSE`.

- absences:

  Logical. Should absence records be kept? Default is `FALSE`.

- basis:

  Character. Vector giving the accepted bases of record. The default
  keeps commonly used occurrence-oriented record types and excludes
  specimen- and unknown-based records.

- establishment:

  Character. Vector giving accepted `degreeOfEstablishment` values. The
  default keeps native and broadly established records, together with
  records lacking that field.

- add_infos:

  Optional character. Vector of additional GBIF occurrence fields to
  append to the default output.

- time_period:

  Numeric. Year range to retain. Default is `c(1000, 3000)`. Records
  with missing years are kept.

- identic_xy:

  Logical. Should records with identical longitude-latitude pairs be
  kept? Default is `FALSE`.

- wConverted_xy:

  Logical. Should records that appear to be incorrectly converted from
  degree-minute notation be kept? Default is `TRUE`. If `FALSE`, an
  approximate version of
  [`CoordinateCleaner::cd_ddmm()`](https://ropensci.github.io/CoordinateCleaner/reference/cd_ddmm.html)
  is applied.

- centroids:

  Logical. Should records located on raster centroids be kept? Default
  is `FALSE`. If `FALSE`,
  [`CoordinateCleaner::cd_round()`](https://ropensci.github.io/CoordinateCleaner/reference/cd_round.html)
  is used.

- ntries:

  Numeric. Number of download attempts before giving up after GBIF or
  `rgbif` errors. Default is `10`.

- error_skip:

  Logical. If `TRUE`, return an empty result when all download attempts
  fail.

- occ_samp:

  Numeric. Maximum number of GBIF occurrences sampled per geographic
  tile. Default is `10000`.

- should_use_occ_download:

  Logical. If `TRUE`, use
  [`rgbif::occ_download()`](https://docs.ropensci.org/rgbif/reference/occ_download.html)
  instead of
  [`rgbif::occ_search()`](https://docs.ropensci.org/rgbif/reference/occ_search.html).
  This requires GBIF credentials. Default is `FALSE`.

- occ_download_user:

  Character. GBIF username. Required if
  `should_use_occ_download = TRUE`.

- occ_download_pwd:

  Character. GBIF password. Required if
  `should_use_occ_download = TRUE`.

- occ_download_email:

  Character. GBIF email address. Required if
  `should_use_occ_download = TRUE`.

- verbose:

  Logical. Should progress messages be printed? Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`CoordinateCleaner::cd_round()`](https://ropensci.github.io/CoordinateCleaner/reference/cd_round.html).

## Value

A `getGBIF` object (a data.frame subclass) containing the requested
occurrence data. The object may also carry the attributes `filter_log`,
which records how many rows were removed at each filtering step, and
`no_xy`, which stores records without coordinates when they are
requested.

Taxonomic harmonization is reflected in the output columns rather than
by replacing all names with a single accepted label. In particular,
`scientificName` stores the record-level GBIF name, whereas
`acceptedScientificName`, `acceptedTaxonKey`, and `taxonomicStatus`
expose the accepted GBIF interpretation used by the query.

## Details

\#' @details The function follows the same taxonomic matching logic used
by the GBIF website. Internally, the input name is first resolved
against the GBIF backbone taxonomy. If the matched name is a synonym,
the download uses the corresponding accepted GBIF taxon key rather than
the synonym key itself.

Querying by accepted taxon key means that occurrence records linked to
infra-specific taxa (subspecies, varieties) under that accepted name are
also retrieved, mirroring GBIF website behaviour. The `acceptedTaxonKey`
and `scientificName` columns in the output therefore reflect the
original GBIF record-level assignment and may refer to a subspecies or
variety rather than the species-level input name. Use
`get_status(children = TRUE)` beforehand to inspect which infra-specific
keys fall under the queried taxon concept.

The output preserves the original GBIF taxonomic fields for each record
rather than rewriting them to a single accepted name. In particular,
`scientificName` stores the record-level GBIF name, while
`acceptedScientificName`, `acceptedTaxonKey`, and `taxonomicStatus`
reflect the harmonized taxon concept used for the query.

Records linked to non-taxonomic backbone entries (e.g. BOLD barcode
sequences) may also appear in the output, since GBIF associates them
with the accepted taxon key. These records typically carry
`basisOfRecord = "MATERIAL_SAMPLE"`. Note however that `MATERIAL_SAMPLE`
is not a reliable proxy for sequence-based records across all taxonomic
groups. Users can cross-check returned `acceptedTaxonKey` values against
`get_status(children = TRUE)` to identify any records linked to
non-backbone entries; see
[`get_status()`](https://8ginette8.github.io/gbif.range/reference/get_status.md)
examples.

If the requested extent contains many records, the query is fragmented
into spatial tiles so that individual API calls remain manageable. Tiles
are refined until each contains fewer than 10,000 records, which is
faster and more robust than relying on a single large request.

Post-download filtering can remove records outside the study extent,
duplicated coordinates, absences, unwanted bases of record, records with
low spatial precision, suspicious coordinate conversions, and
raster-centroid records. The `grain` argument controls both
coordinate-uncertainty filtering and the minimum number of coordinate
decimals that must be reported.

Decimal filtering follows these thresholds:  

- if 110 km \> `grain` \\\ge\\ 11 km, coordinates with at least 1
  decimal are kept  

- if 11 km \> `grain` \\\ge\\ 1100 m, coordinates with at least 2
  decimals are kept  

- if 1100 m \> `grain` \\\ge\\ 110 m, coordinates with at least 3
  decimals are kept  

- if 110 m \> `grain` \\\ge\\ 11 m, coordinates with at least 4 decimals
  are kept  

- if 11 m \> `grain` \\\ge\\ 1.1 m, coordinates with at least 5 decimals
  are kept

## References

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P.,
Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate,
soil, and land cover on plant species distribution in the European Alps.
Ecological monographs, 91(2), e01433. 10.1002/ecm.1433

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to
the global biodiversity information facility API. 10.5281/zenodo.6023735

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C.,
Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized
cleaning of occurrence records from biological collection databases.
Methods in Ecology and Evolution, 10(5), 744-751.
10.1111/2041-210X.13152

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version
1.6-7." (2022). Terra - CRAN

## See also

[`get_status()`](https://8ginette8.github.io/gbif.range/reference/get_status.md)
to inspect the accepted name and synonym mapping returned by GBIF;
`rgbif` for more general GBIF retrieval workflows; and
`CoordinateCleaner` for more extensive occurrence cleaning.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download worldwide observations of Panthera tigris
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Plot
countries <- rnaturalearth::ne_countries(
       type = "countries",
       returnclass = "sv"
)
terra::plot(countries, col = "#bcbddc")
graphics::points(
       obs.pt[, c("decimalLongitude","decimalLatitude")],
       pch = 20,
       col = "#238b4550",
       cex = 4
)

# Download worldwide observations of Ailuropoda melanoleuca
# (with a 100km grain, after 1990 and by keeping duplicates and by
# adding the name of the person who collected the panda records)
obs.am <- get_gbif(
       sp_name = "Ailuropoda melanoleuca",
       grain = 100 ,
       duplicates = TRUE,
       time_period = c(1990,3000),
       add_infos = c("recordedBy","issue")
)

# Plot
terra::plot(countries, col = "#bcbddc")
graphics::points(
       obs.am[,c("decimalLongitude","decimalLatitude")],
       pch = 20,
       col = "#238b4550",
       cex = 4
)

# Download worldwide observations of Bison bison
# (with a 1km grain, after 1990, and keeping raster centroids)
obs.pc <- get_gbif(
       sp_name = "Bison bison",
       grain = 1,
       time_period = c(1990,3000),
       centroids = TRUE
)
} # }
```
