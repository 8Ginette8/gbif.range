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
  shapefiles should be stored.

## Value

`NULL`. Files are downloaded and unpacked for side effects.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download every ecoregion dataset listed in ecoreg_list
get_ecoreg()
} # }
```
