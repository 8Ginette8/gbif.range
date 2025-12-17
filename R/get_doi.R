### =========================================================================
### get_doi
### =========================================================================
#' Get a custom DOI for a GBIF filtered dataset
#'
#' A small user friendly wrapper of \code{derived_dataset()}
#' (\code{rgbif} R package), compatible with one or several
#' \code{get_gbif()} outputs.
#' 
#' @param gbifs Object of class getGBIF or a list of several.
#' @param title Character. Chosen title for the GBIF dataset.
#' @param description Character. Description of the GBIF dataset.
#' @param source_url  Character. Link to where the dataset is stored.
#' @param user Character. User's GBIF username.
#' @param pwd Character. User's GBIF password.
#' @param ... Additonnal parameters for \code{derived_dataset()} called
#' from the \code{rgbif} R package.
#' @details See \code{derived_dataset()} from the \code{rgbif} R package.
#' @return One citable GBIF DOI and associated user's information.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the
#' global biodiversity information facility API. 10.5281/zenodo.6023735
#' @seealso The \code{rgbif package} for additional and more general approaches
#' to get GBIF DOI.
#' @example inst/examples/get_doi_help.R
#' @importFrom rgbif derived_dataset
#' @importFrom methods is
#' @export
get_doi <- function(gbifs = NULL,
				   title = NULL,
				   description = NULL,
				   source_url = "https://example.com/", 
				   user = "",
				   pwd = "",
				   ...) {

  ######################################################
  ### Stop messages
  ######################################################


  # Test input
  if (!(methods::is(gbifs,"getGBIF") || methods::is(gbifs,"list"))) {
  	stop("Given 'gbifs' must be a getGBIF object.")
  }

  # Other checks
  check_character_vector(title, "title")
  check_character_vector(description, "description")
  check_character_vector(source_url, "source_url")
  check_character_vector(user, "user")
  check_character_vector(pwd, "pwd")


  ######################################################
  ### Code
  ######################################################

	# If data.frame, transform in list for code homogenisation
	if (methods::is(gbifs,"getGBIF")) {gbifs <- list(gbifs)}

	# Combine everything
	all.obs <- do.call("rbind", gbifs)

	# Pre-compile a bit
	d.target <- table(all.obs$datasetKey)
	d.summary <- data.frame(
		datasetKey = names(d.target),
		count = as.numeric(d.target)
	)

	# Run the doi function from rgbif
	rgbif::derived_dataset(citation_data = d.summary,
						   title = title,
						   description = description,
						   source_url = source_url,
						   user = user,
						   pwd = pwd,
						   ...)
}