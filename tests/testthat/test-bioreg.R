# Copyright (c) 2025
#
# Unit tests for bioreg helper functions in the package.
# These tests are written to be CI-friendly: they avoid network access
# and heavy external dependencies by mocking behavior where appropriate.
#
# Notes for maintainers:
# - Tests use `mockery::stub()` to temporarily replace functions such as
#   `system.file()` or `terra::vect()` so the test suite doesn't require
#   an installed terra or internet access.
# - Temporary directories (via `tempdir()`) are used to isolate filesystem
#   side effects; tests clean up or write minimal fake files as needed.
# - Do NOT load packages at top-level inside test files; use namespaced
#   calls (e.g., `mockery::stub`) so code is not executed outside of
#   `test_that()` blocks.

base_dir <- tempdir()

update_reference <- FALSE

## Use functions namespaced to avoid executing library() at top-level inside tests

context("bioreg functions")

test_that("get_save_dir returns package extdata when installed or falls back to inst/extdata/downloads", {
  # Case when system.file returns a path (simulate by stubbing system.file)
  # We're replacing `system.file()` inside `get_save_dir()` so the function
  # behaves as if the package were installed (system.file returns a real path).
  # This avoids relying on an actual installed package and keeps the test fast
  # and deterministic.
  mockery::stub(get_save_dir, 'system.file', file.path(base_dir, 'extdata', 'downloads'))
  expect_equal(get_save_dir(NULL), file.path(base_dir, 'extdata', 'downloads'))

  # Case when system.file returns empty string -> fallback to working dir
  # Simulate the package not being installed: system.file returns "".
  # In that case we expect `get_save_dir()` to fall back to a path inside
  # the current working directory (inst/extdata/downloads).
  mockery::stub(get_save_dir, 'system.file', "")
  wd <- getwd()
  expect_equal(get_save_dir(NULL), file.path(wd, 'inst', 'extdata', 'downloads'))
})

test_that("get_bioreg handles non-existing bioregions", {
  # Use tempdir for save_dir
  sd <- file.path(tempdir(), "gbif_test_downloads")
  # Ensure directory cleaned
  if (dir.exists(sd)) unlink(sd, recursive = TRUE)

  # Request a non-existing bioregion should error
  # We expect a clear error message listing available bioregions; the
  # pattern check here is intentionally simple to avoid brittleness.
  expect_error(get_bioreg(bioreg_name = "not_a_bioreg", save_dir = sd),
               "Bioregion not_a_bioreg not found")
})

test_that("check_and_get_bioreg calls get_bioreg when needed", {
  sd <- file.path(tempdir(), "gbif_test_downloads2")
  if (dir.exists(sd)) unlink(sd, recursive = TRUE)

  # Stub get_bioreg to avoid downloads and to record calls
  called <- FALSE
  fake_get <- function(bioreg_name, save_dir) {
    called <<- TRUE
    dir.create(file.path(save_dir, bioreg_name), recursive = TRUE)
    # create a fake shp file to satisfy subsequent listing
    writeLines("fake", file.path(save_dir, bioreg_name, "fake.shp"))
    invisible(NULL)
  }

  # Replace the real downloader with our fake function. The fake creates a
  # directory and writes a minimal fake .shp filename so the rest of the
  # code (which checks for .shp presence) behaves as expected.
  mockery::stub(check_and_get_bioreg, 'get_bioreg', fake_get)

  # Call the function under test and assert that our fake downloader was
  # invoked and that the expected .shp file exists in the save directory.
  check_and_get_bioreg(bioreg_name = "eco_terra", save_dir = sd)
  expect_true(called)
  expect_true(dir.exists(file.path(sd, "eco_terra")))
  expect_true(length(list.files(file.path(sd, "eco_terra"), pattern = "\\.shp$")) == 1)
})

test_that("read_bioreg returns terra vect object when shp exists", {
  sd <- file.path(tempdir(), "gbif_test_downloads3")
  if (dir.exists(sd)) unlink(sd, recursive = TRUE)
  dir.create(file.path(sd, "eco_terra"), recursive = TRUE)
  shpfile <- file.path(sd, "eco_terra", "fake.shp")
  writeLines("fake", shpfile)

  # mock terra::vect to avoid depending on terra
  # Instead of loading terra and reading a real shapefile, stub terra::vect
  # to return a simple S3 object. This keeps the test fast and avoids the
  # need for the terra dependency in CI for this unit test.
  fake_vect <- function(x) { structure(list(path = x), class = "SpatVector") }
  mockery::stub(read_bioreg, 'terra::vect', fake_vect)

  # Call the loader and assert it returns the shape-like object with the
  # expected path information that points to the fake .shp created above.
  res <- read_bioreg(bioreg_name = "eco_terra", save_dir = sd)
  expect_s3_class(res, "SpatVector")
  expect_equal(res$path, shpfile)
})
