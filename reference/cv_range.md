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
# Load available ecoregions
eco_terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir(),
    format = "sf"
)

# First download the worldwide observations of the panda from GBIF
obs_am <- get_gbif(sp_name = "Ailuropoda melanoleuca")
#> |--------------------------------------------|
#> | Total number (all records)    :        290 |
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

# Both calls above depend on remote services and return NULL or an empty
# table if those are unavailable, so guard the rest of the example
if (gbif_have(eco_terra, obs_am)) {

# Build a range map from occurrence points
range_panda <- get_range(
    occ_coord = obs_am,
    ecoreg = eco_terra,
    ecoreg_name = "ECO_NAME",
    format = "sf"
)
am_test <- cv_range(
    range_object = range_panda,
    cv = "block-cv",
    nfolds = 5,
    nblocks = 2
);am_test

}
#> ## Start of computation for species: Ailuropoda melanoleuca ###
#> 11 outlier's from 39 | proportion from total points: 28%
#> ecoregion 1 of 5: Daba Mountains Evergreen Forests
#> ecoregion 2 of 5: Qin Ling Mountains Deciduous Forests
#> ecoregion 3 of 5: Qionglai-Minshan Conifer Forests
#> ecoregion 4 of 5: South China-Vietnam Subtropical Evergreen Forests
#> ecoregion 5 of 5: Southeast Tibet Shrublands And Meadows
#> ## End of computation for species: Ailuropoda melanoleuca ###
#> 6 variables with 5, 5, 5, 5, ... levels: 15625 function evaluations required.
#> ...fold1
#> ...fold2
#> ...fold3
#> ...fold4
#> ...fold5
#> 
#>        TP FA     TA    FP   Precision Sensitivity Specificity         TSS
#> CV1  12.0  0 1210.0 711.0 0.016597510   1.0000000   0.6298803  0.62988027
#> CV2   5.0  1 1179.0 288.0 0.017064846   0.8333333   0.8036810  0.63701431
#> CV3   3.0  2 1366.0 222.0 0.013333333   0.6000000   0.8602015  0.46020151
#> CV4   3.0  2 2130.0 303.0 0.009803922   0.6000000   0.8754624  0.47546239
#> CV5   0.0  5 2548.0  43.0 0.000000000   0.0000000   0.9834041 -0.01659591
#> Mean  4.6  2 1686.6 313.4 0.011359922   0.6066667   0.8305258  0.43719252
# }
```
