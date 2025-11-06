# Tests for small utility functions: optme, make_blocks, make_tiles
# These tests avoid heavy external dependencies by either using pure inputs
# (optme) or stubbing expensive helpers (terra::ext in make_tiles).

context("utility functions")

test_that("optme computes expected penalty for a simple case", {
  # Build a small synthetic example where we can compute the penalty by hand.
  # grps: two groups with initial counts 1 and 2
  grps <- list(1, 2)

  # x: group membership for three items; nms: associated weights/sizes
  x <- c(1, 2, 1)
  nms <- c(10, 20, 30)

  # After aggregation by group: group1 -> 10+30 = 40, group2 -> 20
  # grp initial sums: 1 and 2 -> after addition: 41 and 22
  # tot is total (sum of the final group sizes): 63
  tot <- 63

  # Manually computed penalty: each group compared to tot/2 = 31.5
  # (41 - 31.5)^2 + (22 - 31.5)^2 = 9.5^2 + (-9.5)^2 = 180.5
  expect_equal(optme(x = x, nms = nms, grps = grps, tot = tot), 180.5)
})


test_that("make_blocks handles missing input and simple sampling", {
  # Missing input should error
  expect_error(make_blocks(nfolds = 3, df = data.frame(), npoints = NA),
               "Please supply number of points if no data.frame is supplied")

  # When no df is supplied but npoints is provided, we get a vector of length npoints
  set.seed(42)
  out <- make_blocks(nfolds = 3, df = data.frame(), npoints = 6)
  expect_length(out, 6)
  # All values should be integers in 1:nfolds
  expect_true(all(out %in% 1:3))
})


test_that("make_blocks creates strata for single-column df", {
  # Create a single-column data.frame with enough rows for nblocks
  df <- data.frame(v = seq(1, 8))
  res <- make_blocks(nfolds = 2, df = df, nblocks = 2)

  # Result should have same length as input and contain exactly two strata labels
  expect_length(res, nrow(df))
  expect_true(all(res %in% 1:2))
  expect_equal(length(unique(res)), 2)
})


test_that("make_tiles returns the expected number of tiles when terra::ext is stubbed", {
  # Stub terra::ext to avoid depending on terra. The stub returns a simple
  # extent object (list with xmin/xmax/ymin/ymax) regardless of input. This
  # allows make_tiles to run the partitioning logic without the terra package.
  fake_ext <- function(...) {
    # Return a simple named list that has the fields used by make_tiles
    structure(list(xmin = -180, xmax = 180, ymin = -90, ymax = 90), class = "SpatExtent")
  }

  mockery::stub(make_tiles, 'terra::ext', fake_ext)

  # For Ntiles = 4 we expect 4 polygon strings (sext = FALSE to simplify)
  tiles <- make_tiles(geo = NULL, Ntiles = 4, sext = FALSE)
  expect_true(is.list(tiles) || is.vector(tiles))
  expect_length(tiles, 4)
  # Each tile should be a POLYGON string
  expect_true(all(grepl('^POLYGON', unlist(tiles))))
})
