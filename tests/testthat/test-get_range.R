test_that("get_range returns a range object for offline fixture data", {
  # Work with one species and one enclosing ecoregion to keep the geometry
  # path deterministic while still exercising the clustering and hull code.
  occ <- single_species_occ(species = "Crocuta crocuta", n = 15)
  ecoreg <- make_test_ecoreg(occ)

  set.seed(1)
  range_obj <- get_range(
    occ_coord = occ,
    ecoreg = ecoreg,
    ecoreg_name = "ECO_NAME",
    degrees_outlier = 30,
    clust_pts_outlier = 2,
    buff_width_point = 1,
    buff_incrmt_pts_line = 0.1,
    buff_width_polygon = 1,
    dir_temp = tempdir(),
    raster = FALSE,
    format = "SpatVector",
    verbose = FALSE
  )

  # get_range() returns a reference object whose spatial output should be non-empty.
  expect_true(methods::is(range_obj, "getRange"))
  expect_true(inherits(range_obj$rangeOutput, "SpatVector"))
  expect_gte(nrow(as.data.frame(range_obj$rangeOutput)), 1)
})
