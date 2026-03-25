test_that("split_gbif_by_species writes one tab-delimited file per species", {
  input_file <- occ_fixture_path()
  outdir <- file.path(tempdir(), "gbif_split_test")
  unlink(outdir, recursive = TRUE)

  # A small chunk size forces multiple read-write cycles while staying fast.
  split_summary <- split_gbif_by_species(
    input_file = input_file,
    outdir = outdir,
    chunk_size = 7,
    sep_in = "\t",
    sep_out = "\t",
    overwrite = TRUE,
    verbose = FALSE
  )

  # The fixture contains two species, so two per-species files should be written.
  expect_equal(nrow(split_summary), 2)
  expect_true(all(file.exists(split_summary$species_file)))
  expect_true(all(grepl("^occurrences_speciesKey_", basename(split_summary$species_file))))
  expect_equal(sort(split_summary$n_records), c(20, 20))

  first_file <- utils::read.delim(
    split_summary$species_file[1],
    sep = "\t",
    stringsAsFactors = FALSE
  )

  # Output files should retain the requested GBIF columns without row names.
  expect_setequal(
    names(first_file),
    c(
      "speciesKey",
      "species",
      "scientificName",
      "decimalLongitude",
      "decimalLatitude"
    )
  )
})
