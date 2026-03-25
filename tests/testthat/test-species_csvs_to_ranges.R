test_that("species_csvs_to_ranges builds offline ranges from split occurrence files", {
  input_file <- occ_fixture_path()
  split_dir <- file.path(tempdir(), "gbif_species_split")
  occ_dir <- file.path(tempdir(), "gbif_occ_min")
  range_dir <- file.path(tempdir(), "gbif_ranges")
  unlink(split_dir, recursive = TRUE)
  unlink(occ_dir, recursive = TRUE)
  unlink(range_dir, recursive = TRUE)

  split_summary <- split_gbif_by_species(
    input_file = input_file,
    outdir = split_dir,
    chunk_size = 9,
    sep_in = "\t",
    sep_out = "\t",
    overwrite = TRUE,
    verbose = FALSE
  )

  # One enclosing polygon is enough to exercise the batch range path offline.
  ecoreg <- make_test_ecoreg(load_occ_fixture(), field = "ECO_NAME")

  range_summary <- species_csvs_to_ranges(
    species_dir = split_dir,
    ecoreg = ecoreg,
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
    raster = FALSE,
    format = "SpatVector",
    verbose = FALSE
  )

  expect_equal(nrow(range_summary), nrow(split_summary))
  expect_true(all(file.exists(range_summary$range_file)))
  expect_true(all(file.exists(range_summary$occ_file)))

  saved_obj <- readRDS(range_summary$range_file[1])
  restored_obj <- read_range_rds(range_summary$range_file[1])

  # Plain readRDS() should yield the packed plot-ready wrapper, while the
  # helper restores the underlying terra object for analysis.
  expect_true(inherits(saved_obj$rangeOutput, "gbifPackedSpatVector"))
  expect_true(inherits(restored_obj$rangeOutput, "SpatVector"))
  expect_gte(terra::nrow(restored_obj$rangeOutput), 1)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  expect_no_error(plot(saved_obj$rangeOutput))
})

test_that("species_csvs_to_ranges resolves built-in ecoregion names", {
  input_file <- occ_fixture_path()
  split_dir <- file.path(tempdir(), "gbif_species_split_shortcut")
  range_dir <- file.path(tempdir(), "gbif_ranges_shortcut")
  unlink(split_dir, recursive = TRUE)
  unlink(range_dir, recursive = TRUE)

  split_gbif_by_species(
    input_file = input_file,
    outdir = split_dir,
    chunk_size = 12,
    sep_in = "\t",
    sep_out = "\t",
    overwrite = TRUE,
    verbose = FALSE
  )

  ecoreg <- make_test_ecoreg(load_occ_fixture(), field = "ECO_NAME")
  testthat::local_mocked_bindings(
    read_ecoreg = function(ecoreg_name, format = "SpatVector", save_dir = NULL) {
      ecoreg
    },
    .package = "gbif.range"
  )

  range_summary <- species_csvs_to_ranges(
    species_dir = split_dir,
    ecoreg = "eco_terra",
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
    raster = FALSE,
    format = "SpatVector",
    verbose = FALSE
  )

  expect_equal(nrow(range_summary), 2)
  expect_true(all(file.exists(range_summary$range_file)))
})
