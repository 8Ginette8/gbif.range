# Tests for get_doi: ensure it builds the citation_data correctly and calls rgbif::derived_dataset

context("get_doi")

test_that("get_doi builds citation_data and calls rgbif::derived_dataset with expected args", {
  # Create a small fake gbifs data.frame with datasetKey values
  df <- data.frame(datasetKey = c('A','A','B'), value = 1:3, stringsAsFactors = FALSE)

  # Prepare a fake derived_dataset that captures its arguments
  called <- FALSE
  captured <- NULL
  fake_derived <- function(...) {
    called <<- TRUE
    captured <<- list(...)
    return(list(doi = "10.000/fake", status = "CREATED"))
  }

  # Stub the rgbif::derived_dataset call inside get_doi
  mockery::stub(get_doi, 'rgbif::derived_dataset', fake_derived)

  # Call get_doi with a data.frame (function should convert to list internally)
  res <- get_doi(gbifs = df, title = "T", description = "D", source_url = "S", user = "u", pwd = "p")

  # Check fake derived_dataset was invoked and return value passed through
  expect_true(called)
  expect_type(res, 'list')
  expect_equal(res$doi, "10.000/fake")

  # Inspect captured citation_data: it should be a data.frame with counts per datasetKey
  expect_true('citation_data' %in% names(captured))
  cit <- captured$citation_data
  # Citation data may be provided in different shapes depending on the caller:
  # - a data.frame with columns datasetKey and count
  # - a table/named integer vector where names are dataset keys
  # - a data.frame with counts and rownames set to dataset keys
  if (is.data.frame(cit)) {
    if (all(c('datasetKey','count') %in% names(cit))) {
      counts <- setNames(as.numeric(cit$count), as.character(cit$datasetKey))
    } else if (!is.null(rownames(cit)) && ncol(cit) >= 1 && is.numeric(cit[[1]])) {
      counts <- setNames(as.numeric(cit[[1]]), rownames(cit))
    } else {
      stop('unexpected data.frame structure for citation_data')
    }
  } else if (is.table(cit) || (is.numeric(cit) && !is.null(names(cit)))) {
    counts <- as.numeric(cit)
    names(counts) <- names(cit)
  } else {
    stop('unexpected citation_data type in captured args')
  }

  # Citation data should contain datasetKey 'A' with count 2 and 'B' with count 1
  expect_true(all(c('A','B') %in% names(counts)))
  expect_equal(as.numeric(counts['A']), 2)
  expect_equal(as.numeric(counts['B']), 1)
})
