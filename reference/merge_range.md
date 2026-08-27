# Merge the polygons of a species range map

[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
returns one polygon per occupied ecoregion. Those polygons do not
overlap, so their areas are additive and the ecoregion origin of every
piece stays visible. `merge_range()` dissolves them into a single
polygon, or into groups given by `by`. The dissolve is kept out of
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
because collapsing the ecoregions discards the `ecoreg_name` attribute.
Rasterized output (`format = "SpatRaster"`) never needs this function.

## Usage

``` r
merge_range(x, by = NULL, format = c("SpatVector", "sf"))
```

## Arguments

- x:

  A `getRange` object returned by
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md),
  or the `SpatVector` or `sf` object held in its `rangeOutput` field.

- by:

  Character. Optional name of a field in `x` used to group the polygons:
  one feature is returned per distinct value. Note that on the default
  output of
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
  this is already the case, so `by = ecoreg_name` returns the input
  unchanged. It is useful for coarser groupings, such as a biome field,
  or for objects whose polygons were not built by
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md).
  When `NULL` (default), every polygon is dissolved into a single
  feature.

- format:

  Character. Output format. One of `"SpatVector"` (default) or `"sf"`.

## Value

A `SpatVector` or `sf` object holding the dissolved range. With `by`
supplied, the grouping field is retained along with any other field that
is constant within each group. With `by = NULL` the result carries no
attributes, since no attribute of the original polygons describes the
single merged geometry.

## See also

[`get_range`](https://8ginette8.github.io/gbif.range/reference/get_range.md)()
to build the range map.

## Examples

``` r
# \donttest{
# Load available ecoregions
eco_terra <- read_ecoreg(ecoreg_name = "eco_terra", save_dir = tempdir())

# Worldwide observations of the giant panda from GBIF
obs_am <- get_gbif(sp_name = "Ailuropoda melanoleuca")
#> |--------------------------------------------|
#> | Total number (all records)    :        300 |
#> | Kept records                  :         66 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE
#> 
#> ...GBIF records of Ailuropoda melanoleuca: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> ---------------------------------------------
#>                     step removed remaining
#>          Grain filtering       6        60
#>       Duplicated records      13        47
#>          Absence records       0        47
#>          Basis selection       8        39
#>  Establishment selection       0        39
#>               Time frame       0        39
#>        Identical records       0        39
#>         Raster centroids       0        39
#> 
#> Initial records         : 66
#> Total removed           : 27
#> Final records (XY)      : 39
#> ---------------------------------------------
#> Final records (no XY)   : 0

# Guard against an unavailable remote source
if (gbif_have(eco_terra, obs_am)) {

# Range map: one polygon per occupied ecoregion
rng <- get_range(occ_coord = obs_am, ecoreg = eco_terra,
                 ecoreg_name = "ECO_NAME")
nrow(rng$rangeOutput)
sum(terra::expanse(rng$rangeOutput, unit = "km"))

# Dissolve into a single polygon. The area is unchanged, which shows the
# ecoregion polygons did not overlap.
merged <- merge_range(rng)
nrow(merged)
sum(terra::expanse(merged, unit = "km"))

}
#> ## Start of computation for species: Ailuropoda melanoleuca ###
#> 11 outlier's from 39 | proportion from total points: 28%
#> ecoregion 1 of 5: Daba Mountains Evergreen Forests
#> ecoregion 2 of 5: Qin Ling Mountains Deciduous Forests
#> ecoregion 3 of 5: Qionglai-Minshan Conifer Forests
#> ecoregion 4 of 5: South China-Vietnam Subtropical Evergreen Forests
#> ecoregion 5 of 5: Southeast Tibet Shrublands And Meadows
#> ## End of computation for species: Ailuropoda melanoleuca ###
#> [1] 723407.5
# }
```
