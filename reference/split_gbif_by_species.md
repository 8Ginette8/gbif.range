# Split a Downloaded GBIF Table into One File per Species

Stream a large GBIF export from disk in chunks and write one occurrence
file per species or GBIF taxon key. The function is designed for
multi-species tables that are too large to load fully into memory.

## Usage

``` r
split_gbif_by_species(
  input_file,
  outdir = file.path(tempdir(), "gbif_by_species"),
  chunk_size = 1e+05,
  select_cols = c("speciesKey", "species", "scientificName", "decimalLongitude",
    "decimalLatitude"),
  sep_in = "\t",
  sep_out = "\t",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- input_file:

  Character. Path to a tabular GBIF export already stored on disk.

- outdir:

  Character. Directory where one per-species file will be written.

- chunk_size:

  Integer. Number of rows read at a time. Larger values are usually
  faster, whereas smaller values reduce peak memory use.

- select_cols:

  Character. Vector of columns to keep from the original file. The
  defaults retain the taxon key, species labels, and geographic
  coordinates needed by downstream range workflows.

- sep_in:

  Character. Field separator Used by the input file. GBIF downloads are
  usually tab-delimited.

- sep_out:

  Character. Field separator used for the saved species files.

- overwrite:

  Logical. If `TRUE`, existing batch files created by this function in
  `outdir` are removed before writing new ones.

- verbose:

  Logical. Should progress messages be printed?

## Value

A data frame summarizing the written files, with one row per species key
and the columns `species_key`, `species_name`, `n_records`, and
`species_file`.

## Details

Output file names follow the pattern
`occurrences_speciesKey_<key>_<species>.csv`. The extension is kept as
`.csv` for convenience, even when the file remains tab-delimited.

## See also

[`species_csvs_to_ranges()`](https://8ginette8.github.io/gbif.range/reference/species_csvs_to_ranges.md)
to process the written species files sequentially with
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md).

## Examples

``` r
if (requireNamespace("data.table", quietly = TRUE)) {
  gbif_file <- system.file("extdata", "occ_example_2sps.csv", package = "gbif.range")
  split_dir <- file.path(tempdir(), "gbif_split_help")

  # Remove earlier temporary outputs so the example can be rerun cleanly.
  unlink(split_dir, recursive = TRUE)

  split_summary <- split_gbif_by_species(
    input_file = gbif_file,
    outdir = split_dir,
    chunk_size = 10,
    sep_in = "\t",
    sep_out = "\t",
    overwrite = TRUE,
    verbose = FALSE
  )

  split_summary[, c("species_name", "n_records", "species_file")]
}
#>      species_name n_records
#> 2 Crocuta crocuta        20
#> 1   Hyaena hyaena        20
#>                                                                         species_file
#> 2 /tmp/RtmpXgHW0u/gbif_split_help/occurrences_speciesKey_5218781_Crocuta_crocuta.csv
#> 1   /tmp/RtmpXgHW0u/gbif_split_help/occurrences_speciesKey_5218777_Hyaena_hyaena.csv
```
