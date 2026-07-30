# Read an Ecoregion Layer

Load one of the packaged ecoregion datasets listed in
[`ecoreg_list()`](https://8ginette8.github.io/gbif.range/reference/ecoreg_list.md).
If the requested data are not available locally, they are downloaded
first.

## Usage

``` r
read_ecoreg(ecoreg_name = "eco_terra", save_dir = NULL, format = "SpatVector")
```

## Arguments

- ecoreg_name:

  Character. File name of the ecoregion dataset to load. See
  `ecoreg_list` for available options.

- save_dir:

  Character. Directory where the downloaded files are stored. Defaults
  to [`tempdir()`](https://rdrr.io/r/base/tempfile.html) unless the
  `gbif.range.save_dir` option or `GBIF_RANGE_SAVE_DIR` environment
  variable is set, in which case that location is used instead.

- format:

  `"SpatVector"` (default) or `"sf"`.

## Value

An ecoregion layer as a `SpatVector` or `sf` object, depending on
`format`.

## Details

Four datasets are currently available:

\(1\) `eco_terra` for terrestrial species, based on The Nature
Conservancy version adapted from Olson et al. (2001).

\(2\) `eco_marine` for marine species, based on Spalding et al. (2007,
2012).

\(3\) `eco_hd_marine`, a higher-resolution marine version.

\(4\) `eco_fresh` for freshwater species, based on Abell et al. (2008).

## References

Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand,
H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H.,
Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R.
2001. Terrestrial ecoregions of the world: a new map of life on Earth.
BioScience 51(11):933-938.
[doi:10.1641/0006-3568(2001)051\[0933:TEOTWA\]2.0.CO;2](https://doi.org/10.1641/0006-3568%282001%29051%5B0933%3ATEOTWA%5D2.0.CO%3B2)

The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
Units. GIS layers developed by The Nature Conservancy with multiple
partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986.
Cambridge (UK): The Nature Conservancy.

Spalding, M. D., Fox, H. E., Allen, G. R., Davidson, N., Ferdana, Z. A.,
Finlayson, M., Halpern, B. S., Jorge, M. A., Lombana, A., Lourie, S. A.,
Martin, K. D., McManus, E., Molnar, J., Recchia, C. A., Robertson, J.
(2007). Marine Ecoregions of the World: A Bioregionalization of Coastal
and Shelf Areas. BioScience, 57(7), 573-583.
[doi:10.1641/B570707](https://doi.org/10.1641/B570707)

Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
Pelagic provinces of the world: a biogeographic classification of the
world's surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
[doi:10.1016/j.ocecoaman.2011.12.016](https://doi.org/10.1016/j.ocecoaman.2011.12.016)

The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
of the World. GIS layers developed by The Nature Conservancy with
multiple partners, combined from Spalding et al. (2007) and Spalding et
al. (2012). Cambridge (UK): The Nature Conservancy.

Abell, R., Thieme, M. L., Revenga, C., Bryer, M., Kottelat, M.,
Bogutskaya, N., Coad, B., Mandrak, N., Contreras Balderas, S., Bussing,
W., Stiassny, M. L. J., Skelton, P., Allen, G. R., Unmack, P., Naseka,
A., Ng, R., Sindorf, N., Robertson, J., Armijo, E., Higgins, J. V.,
Heibel, T. J., Wikramanayake, E., Olson, D., Lopez, H. L., Reis, R. E.,
Lundberg, J. G., Sabaj Perez, M. H., Petry, P. (2008). Freshwater
Ecoregions of the World: A New Map of Biogeographic Units for Freshwater
Biodiversity Conservation. BioScience, 58(5), 403-414.
[doi:10.1641/B580507](https://doi.org/10.1641/B580507)

## See also

[`ecoreg_list`](https://8ginette8.github.io/gbif.range/reference/ecoreg_list.md)
for the list of available datasets,
[`get_ecoreg`](https://8ginette8.github.io/gbif.range/reference/get_ecoreg.md)()
to force a fresh download, and
[`check_and_get_ecoreg`](https://8ginette8.github.io/gbif.range/reference/check_and_get_ecoreg.md)()
for the underlying download check.

## Examples

``` r
# \donttest{
shp_eco_terra <- read_ecoreg("eco_terra", save_dir = tempdir())
terra::plot(shp_eco_terra)

# }
```
