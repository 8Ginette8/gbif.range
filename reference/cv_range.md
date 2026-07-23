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
eco.terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir(),
    format = "sf"
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pd <- get_gbif(sp_name = "Ailuropoda melanoleuca")
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
#>          Basis selection      10        37
#>  Establishment selection       0        37
#>               Time frame       0        37
#>        Identical records       0        37
#>         Raster centroids       0        37
#> 
#> Initial records         : 66
#> Total removed           : 29
#> Final records (XY)      : 37
#> ---------------------------------------------
#> Final records (no XY)   : 0

# Build a range map from occurrence points
range.panda <- get_range(
    occ_coord = obs.pd,
    ecoreg = eco.terra,
    clust_pts_outlier = 6,
    ecoreg_name = "ECO_NAME",
    format = "sf"
)
#> ## Start of computation for species:  Ailuropoda melanoleuca  ### 
#> 12 outlier's from 37 | proportion from total points: 32%
#> ecoregion 1  of  4 :  Daba Mountains Evergreen Forests 
#> ecoregion 2  of  4 :  Qin Ling Mountains Deciduous Forests 
#> ecoregion 3  of  4 :  Qionglai-Minshan Conifer Forests 
#> ecoregion 4  of  4 :  Southeast Tibet Shrublands And Meadows 
#> ## End of computation for species:  Ailuropoda melanoleuca  ### 
pd.test <- cv_range(
    range_object = range.panda,
    cv = "block-cv",
    nfolds = 5,
    nblocks = 2
);pd.test
#> 6 variables with 5, 5, 5, 5, ... levels: 15625 function evaluations required.
#> ...fold1
#> ...fold2
#> ...fold3
#> ...fold4
#> ...fold5
#> 
#>       TP  FA     TA     FP   Precision Sensitivity Specificity       TSS
#> CV1  9.0 0.0   66.0  114.0 0.073170732        1.00   0.3666667 0.3666667
#> CV2  5.0 0.0 2071.0  665.0 0.007462687        1.00   0.7569444 0.7569444
#> CV3  5.0 0.0  755.0 1438.0 0.003465003        1.00   0.3442772 0.3442772
#> CV4  3.0 1.0  570.0  640.0 0.004665630        0.75   0.4710744 0.2210744
#> CV5  2.0 2.0 2926.0  755.0 0.002642008        0.50   0.7948927 0.2948927
#> Mean 4.8 0.6 1277.6  722.4 0.018281212        0.85   0.5467711 0.3967711

# }
```
