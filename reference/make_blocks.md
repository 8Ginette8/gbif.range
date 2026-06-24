# Split Data into Approximately Balanced Folds

Create a fold-assignment vector for random or spatially structured
cross-validation.

## Usage

``` r
make_blocks(
  nfolds = 5,
  df = data.frame(),
  nblocks = nfolds * 2,
  npoints = NA,
  pres = numeric()
)
```

## Arguments

- nfolds:

  Numeric. Number of folds to create.

- df:

  Optional `data.frame`. Columns define the structure used to build
  folds. When omitted, random folds are created from `npoints`.

- nblocks:

  Numeric. Number of initial clusters used to build the folds. Must be
  at least `nfolds`.

- npoints:

  Optional numeric. Number of observations to split when `df` is not
  supplied.

- pres:

  Optional binary vector. Used to restrict CLARA clustering to a subset
  of rows and assign the remainder by nearest neighbours.

## Value

An integer vector of length `nrow(df)` or `npoints`, giving the fold
assignment for each observation.

## Details

If `df` has one column, folds are based on quantile bins. If `df` has
two or more columns, folds are based on CLARA clustering. Remaining
clusters are assigned to folds with an optimization step that tries to
balance fold sizes as evenly as possible.

## References

Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O.,
Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species
distribution projections under climate change. Journal of Biogeography,
47(1), 130-142.

## Examples

``` r
if (FALSE) { # \dontrun{
# Downloading worldwide the observations of Panthera tigris
obs.pt <- get_gbif(
    sp_name = "Panthera tigris",
    basis = c("OBSERVATION","HUMAN_OBSERVATION",
        "MACHINE_OBSERVATION","OCCURRENCE")
)

# Create a vector of folds (n = 5) spatially blocked (n = 10)
block.pt <- make_blocks(
    nfolds = 5,
    df = obs.pt[, c("decimalLatitude","decimalLongitude")],
    nblocks = 10
)

# Plot one colour per fold
countries <- rnaturalearth::ne_countries(
    type = "countries",
    returnclass = "sv"
)
countries.focus <- terra::crop(
    countries,
    terra::ext(60,100,0,40)
)
terra::plot(countries.focus, col = "#bcbddc")
graphics::points(
    obs.pt[, c("decimalLongitude","decimalLatitude")],
    pch = 20,
    col = block.pt,
    cex = 1
)
} # }
```
