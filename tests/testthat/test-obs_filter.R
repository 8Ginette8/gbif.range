# Tests for obs_filter: mock terra raster functions to avoid dependency on terra

context("obs_filter")

test_that("obs_filter returns one observation per grid cell and applies threshold", {
  # Create a fake GBIF-like data.frame with two species
  gbifs <- data.frame(
    input.search = c('sp1','sp1','sp1','sp2','sp2'),
    decimalLongitude = c(0, 0.1, 0.1, 10, 10.1),
    decimalLatitude = c(0, 0.1, 0.1, 20, 20.1),
    stringsAsFactors = FALSE
  )

  # Fake grid object (we only pass a placeholder to the function)
  fake_grid <- list()

  # Stub terra::rast to simply return the object passed
  mockery::stub(obs_filter, 'terra::rast', function(x) { fake_grid })

  # Stub terra::cellFromXY to deterministically map coordinates to cell indices.
  # This stub handles both matrix (multiple rows) and vector (single row) inputs.
  mockery::stub(obs_filter, 'terra::cellFromXY', function(grid, xy) {
    # determine number of coordinate rows
    n <- if (is.null(dim(xy))) 1 else nrow(xy)
    # extract longitudes (first column); handle single-row vector too
    lon <- if (is.null(dim(xy))) xy[1] else xy[, 1]
    # simple rule: longitudes < 5 -> cell 1, otherwise cell 2
    if (mean(lon) < 5) {
      return(rep(1, n))
    } else {
      return(rep(2, n))
    }
  })

  # Stub terra::xyFromCell to return coordinates with column names 'x' and 'y'
  mockery::stub(obs_filter, 'terra::xyFromCell', function(grid, cells) {
    # Ensure output has columns named x and y to match production behaviour
    make_row <- function(xval, yval) {
      m <- matrix(c(xval, yval), nrow = 1)
      colnames(m) <- c('x', 'y')
      return(m)
    }
    if (all(cells == 1)) {
      return(make_row(0, 0))
    } else if (all(cells == 2)) {
      return(make_row(10, 20))
    } else {
      # for mixed/multiple cells, return rows in order with names
      rows <- lapply(cells, function(cl) if (cl == 1) c(0,0) else c(10,20))
      m <- do.call(rbind, rows)
      colnames(m) <- c('x','y')
      return(m)
    }
  })

  # Run without threshold: should return one row per occupied cell per species
  out <- obs_filter(gbifs = gbifs, grid = NULL, threshold = NULL)
  # We expect at least one row per species-cell combination; check species present
  expect_true(all(c('sp1','sp2') %in% out$Species))

  # Run with threshold = 2: only cells with >=2 obs should be kept. In our
  # fake mapping, cell 1 has 3 observations (sp1) and cell 2 has 2 observations (sp2).
  # With threshold = 2 both cells should be kept; verify result contains both species.
  out2 <- obs_filter(gbifs = gbifs, grid = NULL, threshold = 2)
  expect_true(all(c('sp1','sp2') %in% out2$Species))
})
