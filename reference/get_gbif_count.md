# Count GBIF Occurrences Before Downloading Data

Query GBIF and return the number of available records for a taxon using
the same name-matching logic as
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md),
but without downloading the occurrence table.

## Usage

``` r
get_gbif_count(
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
  spatial_issue = FALSE
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

  Logical. If `TRUE` (default), count only records with coordinates. If
  `FALSE`, count only records without coordinates. If `NULL`, count all
  records.

- spatial_issue:

  Logical. If `FALSE` (default), count only records without geospatial
  issues. If `TRUE`, count only records with geospatial issues. If
  `NULL`, count all records.

## Value

A numeric vector of length two giving the total number of GBIF records
and the number retained by the requested filters. A formatted summary is
also printed to the console.

## Details

The function mirrors the taxonomic matching strategy used by
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md),
then reports both the total number of GBIF records and the number
retained after applying the chosen filters.

## References

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to
the global biodiversity information facility API. 10.5281/zenodo.6023735

## Examples

``` r
if (FALSE) { # \dontrun{
# Get number of observations with default filters
obs.pt <- get_gbif_count(
  sp_name = "Ailuropoda melanoleuca",
  has_xy = TRUE,
  spatial_issue = FALSE,
  geo = NULL
)

# Get the total number of observations
obs.pt <- get_gbif_count(
  sp_name = "Ailuropoda melanoleuca",
  has_xy = NULL,
  spatial_issue = NULL,
  geo = NULL
)

# Example of setting global 'geo' (all records are still kept)
obs.pt <- get_gbif_count(
  sp_name = "Ailuropoda melanoleuca",
  has_xy = NULL,
  spatial_issue = NULL,
  geo = terra::ext()
)

# Example of fuzzy matching when search is set to FALSE
obs.pt <- get_gbif_count(
  sp_name = "Ailuropoda melanolca",
  search = FALSE
)

# Example on the European Alps
shp.lonlat <- terra::vect(
    paste0(
        system.file(package = "gbif.range"),
        "/extdata/shp_lonlat.shp"
    )
)
obs.pt <- get_gbif_count(
  sp_name = "Arctostaphylos alpinus",
  has_xy = TRUE,
  spatial_issue = FALSE,
  geo = shp.lonlat
)

} # }
```
