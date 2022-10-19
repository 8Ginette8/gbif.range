### =========================================================================
### wsl.gbif
### =========================================================================
#' Get a custom DOI for a GBIF filtered dataset
#'
#' A small user friendly wrapper of the derived_dataset() function of the
#' rgbif R package, compatible with one or several wsl_gbif() outputs.
#' 
#' 
#' @param wsl_gbif List. List of one or several wsl_gbif outputs.
#' @param title The title for your derived dataset.
#' @param descritpion A description of the dataset.
#' @param source_url  A link to where the dataset is stored.
#' @param user Your GBIF username.
#' @param pwd Your GBIF password.
#' @param ... Additonnal parameters for derived_dataset() in rgbif.
#' R package.
#' @details see derived_dataset() function from the rgbif R package
#' @return One citable DOI and its information.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches to get GBIF DOI
#' @examples
#' 
#' # Necessary libraries
#' #library(rgbif)
#' 
#' # Downloading worldwide the observations of Panthera tigris and Ailuropoda melanoleuca
#' test1 = wsl_gbif("Panthera tigris")
#' test32 = wsl_gbif("Ailuropoda melanoleuca")
#' 
#' # Just an example on how to retrieve the DOI for only one dataset
#' 
#' 
#' # Just an example on how to retrieve the DOI for only one dataset
#' 
#' 
#' d.target = table(test1$datasetKey)
#' d.summary = data.frame(datasetKey = names(d.target),count = as.numeric(d.target))
#' rgbif::derived_dataset(d.summary,"GBIF_test",
#'     "Filetred and cleaned based on CoordinateCleaner",source_url="https://example.com/",
#'     user="your_gbif_user",pwd="your_gbif_password")
#' 
#' @export
wsl_doi = function(wsl_gbif = list(),
				   title = NULL,
				   description = NULL,
				   source_url = "https://example.com/", 
				   usr = "",
				   pwd = "",
				   ...) {

}