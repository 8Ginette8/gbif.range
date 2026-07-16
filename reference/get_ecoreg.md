# Download Ecoregion Layers

Download one or more ecoregion datasets listed in
[`ecoreg_list()`](https://8ginette8.github.io/gbif.range/reference/ecoreg_list.md)
and unpack them into a local directory.

## Usage

``` r
get_ecoreg(ecoreg_name = "all", save_dir = NULL)
```

## Arguments

- ecoreg_name:

  Character. Use `"all"` to download every dataset, or supply a single
  file name listed in `ecoreg_list`.

- save_dir:

  Character. Directory where the downloaded zip files and extracted
  shapefiles should be stored. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) unless the
  `gbif.range.save_dir` option or `GBIF_RANGE_SAVE_DIR` environment
  variable is set, in which case that location is used instead.

## Value

`NULL`. Files are downloaded and unpacked for side effects.

## See also

[`ecoreg_list`](https://8ginette8.github.io/gbif.range/reference/ecoreg_list.md)
for the list of available datasets, and
[`read_ecoreg`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)()
to download and load a layer in one step.

## Examples

``` r
# \donttest{
# Download every ecoregion dataset listed in ecoreg_list
get_ecoreg(save_dir = tempdir())
#> Preparing to download all listed ecoregions: eco_terra ; eco_fresh ; eco_marine ; eco_hd_marine 
#>  data files to: /tmp/RtmpXgHW0u
#> Downloaded: eco_terra
#> Description: Terrestrial Ecoregions of the World
#> Unzipped: eco_terra 
#>  saved to: /tmp/RtmpXgHW0u/eco_terra 
#>  removed:  eco_terra.zip
#> Downloaded: eco_fresh
#> Description: Freshwater Ecoregions of the World
#> Unzipped: eco_fresh 
#>  saved to: /tmp/RtmpXgHW0u/eco_fresh 
#>  removed:  eco_fresh.zip
#> Downloaded: eco_marine
#> Description: Marine Ecoregions of the World
#> Unzipped: eco_marine 
#>  saved to: /tmp/RtmpXgHW0u/eco_marine 
#>  removed:  eco_marine.zip
#> Downloaded: eco_hd_marine
#> Description: High Definition Marine Ecoregions and Pelagic Provinces of the World (2007, 2012)
#> Unzipped: eco_hd_marine 
#>  saved to: /tmp/RtmpXgHW0u/eco_hd_marine 
#>  removed:  eco_hd_marine.zip
# }
```
