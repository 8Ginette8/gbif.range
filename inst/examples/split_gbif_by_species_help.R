\dontrun{
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
}
