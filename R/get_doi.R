### =========================================================================
### wsl.gbif
### =========================================================================
#' Get a custom DOI for a GBIF filtered dataset
#'
#' A small user friendly wrapper of the derived_dataset() function of the
#' rgbif R package, compatible with one or several wsl_gbif() outputs.
#' 
#' 
#' @param wsl_gbif data.frame or list. One wsl_gbif() output or a list of several.
#' @param title The title for your derived dataset.
#' @param descritpion A description of the dataset.
#' @param source_url  A link to where the dataset is stored.
#' @param user Your GBIF username.
#' @param pwd Your GBIF password.
#' @param ... Additonnal parameters for derived_dataset() in rgbif.
#' @details see derived_dataset() function from the rgbif R package
#' @return One citable DOI and its information.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches to get GBIF DOI
#' @examples
#' 
#' # Downloading worldwide the observations of Panthera tigris and Ailuropoda melanoleuca
#' obs.pt = wsl_gbif("Panthera tigris")
#' obs.am = wsl_gbif("Ailuropoda melanoleuca")
#' 
#' # Just an example on how to retrieve the DOI for only one wsl_gbif() output
#' wsl_doi(obs.pt,title="GBIF_test1",description="A small example 1",
#'     source_url="https://example.com/",user="",pwd="") # Use your own GBIF credentials here
#' 
#' # Just an example on how to retrieve the DOI for several wsl_gbif() outputs
#' wsl_doi(list(obs.pt,obs.am),title="GBIF_test2",description="A small example 2",
#'     source_url="https://example.com/",user="",pwd="") # Use your own GBIF credentials here
#' 
#' @export
wsl_doi = function(wsl_gbif = NULL,
				   title = NULL,
				   description = NULL,
				   source_url = "https://example.com/", 
				   user = "",
				   pwd = "",
				   ...) {

	# If data.frame, transform in list for code homogenisation
	if (class(wsl_gbif)%in%"data.frame") {wsl_gbif=list(wsl_gbif)}

	# Combine everything
	all.obs = do.call("rbind",wsl_gbif)

	# Pre-compile a bit
	d.target = table(all.obs$datasetKey)
	d.summary = data.frame(datasetKey=names(d.target),count=as.numeric(d.target))

	# Run the doi function from rgbif
	rgbif::derived_dataset(citation_data=d.summary,title=title,
		description=description,source_url=source_url,user=user,pwd=pwd,...)
}