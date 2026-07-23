# Check for a Local Ecoregion Layer and Download It if Needed

Check whether a target ecoregion directory exists and contains at least
one shapefile. If not, download the dataset with
[`get_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/get_ecoreg.md).

## Usage

``` r
check_and_get_ecoreg(ecoreg_name = "eco_terra", save_dir = NULL)
```

## Arguments

- ecoreg_name:

  Character. File name of the ecoregion dataset to check. See
  `ecoreg_list` for valid values.

- save_dir:

  Character. Directory where downloaded files should be stored. Defaults
  to [`tempdir()`](https://rdrr.io/r/base/tempfile.html) unless the
  `gbif.range.save_dir` option or `GBIF_RANGE_SAVE_DIR` environment
  variable is set, in which case that location is used instead.

## Value

`NULL`. The function downloads data only when required.

## See also

[`get_ecoreg`](https://8ginette8.github.io/gbif.range/reference/get_ecoreg.md)()
to force a download, and
[`read_ecoreg`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)()
which calls this function internally before loading the layer.

## Examples

``` r
# \donttest{
check_and_get_ecoreg("eco_terra", save_dir = tempdir())
#> eco_terra directory does not exist or contains no .shp files. 
#>             
#>  [/tmp/RtmpNH8Rzr/eco_terra/eco_terra] will be created and data will be downloaded.
#>             
#>  Downloading data...
#> Preparing to download ecoregion eco_terra 
#>  data file to: /tmp/RtmpNH8Rzr
#> Downloaded: eco_terra
#> Description: Terrestrial Ecoregions of the World
#> Unzipped: eco_terra 
#>  saved to: /tmp/RtmpNH8Rzr/eco_terra 
#>  removed:  eco_terra.zip
# }
```
