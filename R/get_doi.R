### =========================================================================
### get_doi
### =========================================================================
#' Get a custom DOI for a GBIF filtered dataset
#'
#' A small user friendly wrapper of the derived_dataset() function of the
#' rgbif R package, compatible with one or several get_gbif() outputs.
#' 
#' 
#' @param gbifs data.frame or list. One get_gbif() output or a list of several.
#' @param title Title for your derived dataset.
#' @param description Description of the dataset.
#' @param source_url  Link to where the dataset is stored.
#' @param user Your GBIF username.
#' @param pwd Your GBIF password.
#' @param ... Additonnal parameters for derived_dataset() in rgbif.
#' @details see derived_dataset() function from the rgbif R package.
#' @return One citable DOI and its information.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches to get GBIF DOI.
#' @example inst/examples/get_doi_help.R
#' @importFrom rgbif derived_dataset
#' @export
get_doi = function(gbifs = NULL,
				   title = NULL,
				   description = NULL,
				   source_url = "https://example.com/", 
				   user = "",
				   pwd = "",
				   ...) {

	# If data.frame, transform in list for code homogenisation
	if (class(gbifs)%in%"data.frame") {gbifs=list(gbifs)}

	# Combine everything
	all_obs = do.call("rbind",gbifs)

	# Pre-compile a bit
	d_target = table(all_obs$datasetKey)
	d_summary = data.frame(datasetKey=names(d_target),count=as.numeric(d_target))

	# Run the doi function from rgbif
	rgbif::derived_dataset(citation_data=d_summary,title=title,
		description=description,source_url=source_url,user=user,pwd=pwd,...)
}