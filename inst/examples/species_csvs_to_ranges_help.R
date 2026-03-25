\dontrun{
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

  range_summary[, c("species_name", "n_points", "range_file")]

  # The same batch call can also resolve the built-in layer internally:
  # species_csvs_to_ranges(species_dir = split_dir, ecoreg = "eco_terra", ...)
}
}
