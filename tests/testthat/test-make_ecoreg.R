test_that("make_ecoreg creates offline polygon ecoregions", {
  # Two small rasters are enough to exercise the CLARA-based partitioning code.
  r1 <- terra::rast(
    ncols = 4,
    nrows = 4,
    xmin = 0,
    xmax = 4,
    ymin = 0,
    ymax = 4,
    crs = "EPSG:4326"
  )
  r2 <- r1
  terra::values(r1) <- rep(c(1, 2, 3, 4), each = 4)
  terra::values(r2) <- rep(c(4, 3, 2, 1), times = 4)
  env <- c(r1, r2)

  set.seed(1)
  ecoreg <- make_ecoreg(env = env, nclass = 2, raster = FALSE)

  # Polygon output should include the fields created by make_ecoreg().
  expect_true(inherits(ecoreg, "SpatVector"))
  expect_true(all(c("CLARA", "EcoRegion") %in% names(ecoreg)))
  expect_gte(nrow(as.data.frame(ecoreg)), 1)
})
