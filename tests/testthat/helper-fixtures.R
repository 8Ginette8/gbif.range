occ_fixture_path <- function() {
  testthat::test_path("fixtures", "occ_example_2sps.csv")
}

load_occ_fixture <- function() {
  utils::read.delim(
    occ_fixture_path(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

make_test_gbif <- function(x = load_occ_fixture()) {
  x$input_search <- x$species
  gbif.range:::getGBIF(x)
}

single_species_occ <- function(species = "Crocuta crocuta", n = 15) {
  x <- load_occ_fixture()
  x <- x[x$species %in% species, , drop = FALSE]
  if (!is.null(n)) {
    x <- utils::head(x, n)
  }
  x$input_search <- species
  gbif.range:::getGBIF(x)
}

make_test_ecoreg <- function(occ, field = "ECO_NAME", label = "fixture_ecoreg") {
  xmin <- min(occ$decimalLongitude) - 5
  xmax <- max(occ$decimalLongitude) + 5
  ymin <- min(occ$decimalLatitude) - 5
  ymax <- max(occ$decimalLatitude) + 5

  polygon <- sf::st_polygon(list(matrix(
    c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ),
    ncol = 2,
    byrow = TRUE
  )))

  sf::st_sf(
    setNames(data.frame(label, stringsAsFactors = FALSE), field),
    geometry = sf::st_sfc(polygon, crs = 4326)
  )
}
