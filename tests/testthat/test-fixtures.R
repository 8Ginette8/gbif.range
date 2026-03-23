test_that("occurrence fixture can be wrapped as a getGBIF object", {
  occ <- make_test_gbif()

  expect_s3_class(occ, "getGBIF")
  expect_s3_class(occ, "data.frame")
  expect_true(file.exists(occ_fixture_path()))
  expect_true(all(
    c(
      "speciesKey",
      "species",
      "scientificName",
      "decimalLongitude",
      "decimalLatitude",
      "input_search"
    ) %in% names(occ)
  ))
  expect_equal(length(unique(occ$input_search)), 2)
})
