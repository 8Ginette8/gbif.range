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

  Character. Directory where downloaded files should be stored.

## Value

`NULL`. The function downloads data only when required.

## Examples

``` r
if (FALSE) { # \dontrun{
check_and_get_ecoreg("eco_terra")
} # }
```
