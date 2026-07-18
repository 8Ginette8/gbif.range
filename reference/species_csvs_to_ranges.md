# Build Range Maps from Per-Species GBIF Files Saved on Disk

Read one occurrence file per species, prepare the minimal coordinate
table needed by
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md),
and save one range output per species.

## Usage

``` r
species_csvs_to_ranges(
  species_dir,
  ecoreg,
  ecoreg_name = NULL,
  outdir = file.path(tempdir(), "gbif_ranges"),
  occ_outdir = NULL,
  occ_save_as = c("none", "tsv", "rds"),
  range_save_as = c("rds", "gpkg", "tif"),
  deduplicate = TRUE,
  sep_in = "\t",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- species_dir:

  Character. Directory containing files created by
  [`split_gbif_by_species()`](https://8ginette8.github.io/gbif.range/reference/split_gbif_by_species.md).

- ecoreg:

  Spatial ecoregion layer in WGS84. Accepted classes are
  `SpatialPolygonsDataFrame`, `SpatVector`, and `sf`. This can be a
  downloaded layer from
  [`read_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)
  or a custom layer created by
  [`make_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md).

- ecoreg_name:

  Character. String naming the field in `ecoreg` that defines ecoregion
  categories. For `eco_terra`, for example, `"ECO_NAME"` is the most
  detailed level. When `ecoreg` comes from
  [`make_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md),
  `"EcoRegion"` is used automatically.

- outdir:

  Character. Directory where range files will be saved.

- occ_outdir:

  Optional character. Directory where the minimal occurrence tables
  passed to
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
  will also be saved.

- occ_save_as:

  `"none"`, `"tsv"`, or `"rds"`. File format for minimal occurrence
  tables.

- range_save_as:

  `"rds"`, `"gpkg"`, or `"tif"`. File format for range outputs.

- deduplicate:

  Logical. Should identical longitude-latitude pairs be collapsed before
  range inference?

- sep_in:

  Character. Field separator used by the per-species occurrence files.

- overwrite:

  Logical. If `TRUE`, existing outputs in `outdir`/`occ_outdir` are
  replaced.

- verbose:

  Logical. Should progress messages be printed?

- ...:

  Additional arguments passed to
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md).

## Value

A data frame summarizing the processed files, with one row per range and
the columns `species_key`, `species_name`, `n_points`, `occ_file`, and
`range_file`.

## Details

The function is intended to work seamlessly with
[`split_gbif_by_species()`](https://8ginette8.github.io/gbif.range/reference/split_gbif_by_species.md).
It reads each species file sequentially, adds the single-species
`input_search` column expected by
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md),
and keeps only the coordinate columns needed for range construction.

## See also

[`split_gbif_by_species()`](https://8ginette8.github.io/gbif.range/reference/split_gbif_by_species.md)
and
[`read_range_rds()`](https://8ginette8.github.io/gbif.range/reference/read_range_rds.md).

## Examples

``` r
# \donttest{
if (requireNamespace("data.table", quietly = TRUE)) {
  gbif_file <- system.file("extdata", "occ_example_2sps.csv", package = "gbif.range")
  split_dir <- file.path(tempdir(), "gbif_species_help")
  occ_dir <- file.path(tempdir(), "gbif_occ_help")
  range_dir <- file.path(tempdir(), "gbif_range_help")

  unlink(split_dir, recursive = TRUE)
  unlink(occ_dir, recursive = TRUE)
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

  # Crop the packaged terrestrial ecoregions to the extent of the example
  # occurrences so the help example stays lightweight.
  occ_example <- utils::read.delim(gbif_file, sep = "\t", stringsAsFactors = FALSE)
  eco_terra <- read_ecoreg(
    "eco_terra",
    save_dir = tempdir()
    )
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
    occ_outdir = occ_dir,
    occ_save_as = "tsv",
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

  range_summary[, c("species_name", "n_points", "range_file")]

  # The same batch call can also resolve the built-in layer internally:
  range_summary_builtin <- species_csvs_to_ranges(
    species_dir = split_dir,
    ecoreg = "eco_terra",
    ecoreg_name = "ECO_NAME",
    outdir = range_dir,
    overwrite = TRUE
  )

  range_summary_builtin[, c("species_name", "n_points", "range_file")]
}
#> ## Start of computation for species:  Hyaena hyaena hyaena (Linnaeus, 1758)  ### 
#> 4 outlier's from 20 | proportion from total points: 20%
#> ecoregion 1  of  3 :  Arabian Desert And East Sahero-Arabian Xeric Shrublands 
#> ecoregion 2  of  3 :  Eastern Mediterranean Conifer-Sclerophyllous-Broadleaf Forests 
#> ecoregion 3  of  3 :  Mesopotamian Shrub Desert 
#> ## End of computation for species:  Hyaena hyaena hyaena (Linnaeus, 1758)  ### 
#> Saved range for Hyaena hyaena hyaena (Linnaeus, 1758) to range_speciesKey_5218777_Hyaena_hyaena_hyaena_Linnaeus_1758.rds.
#> ## Start of computation for species:  Crocuta crocuta (Erxleben, 1777)  ### 
#> 6 outlier's from 20 | proportion from total points: 30%
#> ecoregion 1  of  5 :  Drakensberg Montane Grasslands, Woodlands And Forests 
#> ecoregion 2  of  5 :  Northern Acacia-Commiphora Bushlands And Thickets 
#> ecoregion 3  of  5 :  Serengeti Volcanic Grasslands 
#> ecoregion 4  of  5 :  Southern Acacia-Commiphora Bushlands And Thickets 
#> ecoregion 5  of  5 :  Zambezian And Mopane Woodlands 
#> ## End of computation for species:  Crocuta crocuta (Erxleben, 1777)  ### 
#> Saved range for Crocuta crocuta (Erxleben, 1777) to range_speciesKey_5218781_Crocuta_crocuta_Erxleben_1777.rds.
#>                            species_name n_points
#> 1 Hyaena hyaena hyaena (Linnaeus, 1758)       20
#> 2      Crocuta crocuta (Erxleben, 1777)       20
#>                                                                                        range_file
#> 1 /tmp/Rtmpc3Ymaq/gbif_range_help/range_speciesKey_5218777_Hyaena_hyaena_hyaena_Linnaeus_1758.rds
#> 2      /tmp/Rtmpc3Ymaq/gbif_range_help/range_speciesKey_5218781_Crocuta_crocuta_Erxleben_1777.rds

# }
```
