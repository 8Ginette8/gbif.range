% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxonomy.R
\name{get_taxonomy}
\alias{get_taxonomy}
\title{Retrieve from GBIF all scientific names of a specific Taxa}
\usage{
get_taxonomy(sp_name = NULL, conf_match = 90, all = FALSE)
}
\arguments{
\item{sp_name}{Character. Species name from which the user wants to retrieve all existing GBIF names.}

\item{conf_match}{Numeric. From 0 to 100. Determine the confidence
threshold of match of 'sp_name' with the GBIF backbone taxonomy. Default is 90.}

\item{all}{Logical. Default is FALSE. Should all species names be retrieved or only
the accepted name and its synonyms?}
}
\value{
Data.frame with three columns: (1) GBIF taxonomic key, (2) scientificName and
(3) Backbone Taxonomy Status.
}
\description{
Generates, based on a given species name, a list of all its scientific names
(accepted, synonyms) found in the GBIF backbone taxonomy and used to download the data in
get_gbif(). Children and related doubtful names not used to download the data may also be extracted.
}
\examples{
get_taxonomy("Cypripedium calceolus",all=FALSE)
get_taxonomy("Cypripedium calceolus",all=TRUE)

}
\references{
Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
information facility API. 10.5281/zenodo.6023735
}
\seealso{
The rgbif package for additional and more general approaches on how to retrieve
scientific names from the GBIF backbone taxonomy.
}
