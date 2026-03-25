occ_fixture_path <- function() {
  # Use the packaged extdata fixture as the single source of truth.
  inst_path <- system.file("extdata", "occ_example_2sps.csv", package = "gbif.range")
  if (nzchar(inst_path)) {
    return(inst_path)
  }

  # During source-tree tests the package may not yet be installed, so fall back
  # to the repository copy under inst/extdata.
  testthat::test_path("..", "..", "inst", "extdata", "occ_example_2sps.csv")
}

load_occ_fixture <- function() {
  # Preserve original column names because downstream package code expects them,
  # and read the GBIF-style file exactly as a tab-delimited export.
  utils::read.delim(
    occ_fixture_path(),
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

make_test_gbif <- function(x = load_occ_fixture()) {
  # Mimic the structure returned by get_gbif() without making web requests.
  x$input_search <- x$species
  gbif.range:::getGBIF(x)
}

single_species_occ <- function(species = "Crocuta crocuta", n = 15) {
  # Use a small deterministic subset so get_range() stays fast and reproducible.
  x <- load_occ_fixture()
  x <- x[x$species %in% species, , drop = FALSE]
  if (!is.null(n)) {
    x <- utils::head(x, n)
  }
  x$input_search <- species
  gbif.range:::getGBIF(x)
}

make_test_ecoreg <- function(occ, field = "ECO_NAME", label = "fixture_ecoreg") {
  # Build one simple polygon around the occurrences so range tests stay offline.
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
