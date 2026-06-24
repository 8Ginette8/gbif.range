# Read a Range File Saved by `species_csvs_to_ranges()`

Convenience wrapper around
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) for range files saved
with `range_save_as = "rds"`.

## Usage

``` r
read_range_rds(file)
```

## Arguments

- file:

  Path to an `.rds` range file created by
  [`species_csvs_to_ranges()`](https://8ginette8.github.io/gbif.range/reference/species_csvs_to_ranges.md).

## Value

A list with `init.args` and `rangeOutput`, matching the structure
written by the batch workflow.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("data.table", quietly = TRUE)) {
  gbif_file <- system.file("extdata", "occ_example_2sps.csv", package = "gbif.range")
  split_dir <- file.path(tempdir(), "gbif_read_range_split")
  range_dir <- file.path(tempdir(), "gbif_read_range_out")

  unlink(split_dir, recursive = TRUE)
  unlink(range_dir, recursive = TRUE)

  split_gbif_by_species(
    input_file = gbif_file,
    outdir = split_dir,
    chunk_size = 10,
    sep_in = "\t",
    sep_out = "\t",
    overwrite = TRUE,
    verbose = FALSE
  )

  occ_example <- utils::read.delim(gbif_file, sep = "\t", stringsAsFactors = FALSE)
  eco_terra <- read_ecoreg("eco_terra")
  eco_crop <- terra::crop(
    eco_terra,
    terra::ext(
      min(occ_example$decimalLongitude) - 5,
      max(occ_example$decimalLongitude) + 5,
      min(occ_example$decimalLatitude) - 5,
      max(occ_example$decimalLatitude) + 5
    )
  )

  range_summary <- species_csvs_to_ranges(
    species_dir = split_dir,
    ecoreg = eco_crop,
    ecoreg_name = "ECO_NAME",
    outdir = range_dir,
    range_save_as = "rds",
    sep_in = "\t",
    overwrite = TRUE,
    degrees_outlier = 30,
    clust_pts_outlier = 2,
    buff_width_point = 1,
    buff_incrmt_pts_line = 0.1,
    buff_width_polygon = 1,
    format = "SpatVector",
    verbose = FALSE
  )

  rg <- read_range_rds(range_summary$range_file[1])
  terra::plot(rg$rangeOutput)
}
} # }
```
