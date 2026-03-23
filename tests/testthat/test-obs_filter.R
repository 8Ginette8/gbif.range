test_that("obs_filter keeps at most one record per species and grid cell", {
  occ <- make_test_gbif()
  grid <- terra::rast(
    xmin = min(occ$decimalLongitude) - 1,
    xmax = max(occ$decimalLongitude) + 1,
    ymin = min(occ$decimalLatitude) - 1,
    ymax = max(occ$decimalLatitude) + 1,
    resolution = 10,
    crs = "EPSG:4326"
  )

  filtered <- obs_filter(occ, grid)
  cell_id <- terra::cellFromXY(grid, as.matrix(filtered[, c("x", "y")]))

  expect_true(all(c("Species", "x", "y") %in% names(filtered)))
  expect_setequal(unique(filtered$Species), unique(occ$input_search))
  expect_equal(anyDuplicated(paste(filtered$Species, cell_id)), 0L)
})

test_that("obs_filter threshold removes sparsely occupied cells", {
  occ_raw <- load_occ_fixture()[1:4, , drop = FALSE]
  occ_raw$decimalLongitude <- c(10.1, 10.2, 10.3, 30.1)
  occ_raw$decimalLatitude <- c(-5.1, -5.2, -5.3, 5.1)
  occ_raw$input_search <- "Crocuta crocuta"
  occ <- gbif.range:::getGBIF(occ_raw)

  grid <- terra::rast(
    xmin = 0,
    xmax = 40,
    ymin = -10,
    ymax = 10,
    resolution = 10,
    crs = "EPSG:4326"
  )

  filtered_all <- obs_filter(occ, grid)
  filtered_threshold <- obs_filter(occ, grid, threshold = 2)

  expect_equal(nrow(filtered_all), 2)
  expect_equal(nrow(filtered_threshold), 1)
  expect_equal(filtered_threshold$Species, "Crocuta crocuta")
})
