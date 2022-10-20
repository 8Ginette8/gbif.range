### ==================================================================
### wsl.taXnames
### ==================================================================
#' Retrieve from GBIF all scientific names of a specific Taxa
#'
#' Generates, based on a given species name, a list of all its scientific names
#' (accepted, synonyms) found in the GBIF backbone taxonomy to download the data.
#' Children and related doubtful names not use to download the data may also be extracted.
#'
#' @param sp_name Character. Species name from which the user wants to retrieve all existing GBIF names
#' @param conf_match Numeric. From 0 to 100. Determine the confidence
#' threshold of match of 'sp_name' with the GBIF backbone taxonomy. Default is 90.
#' @param all Logical. Default is FALSE. Should all species names be retrieved or only
#' the accepted name and its synonyms?
#' 
#' @return A data.frame with two columns: (1) Names and (2) Backbone Taxonomy Status
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches on how to retrieve
#' scientific names from the GBIF backbone taxonomy.
#' @examples
#' wsl_taxonomy("Cypripedium calceolus",all=FALSE)
#' wsl_taxonomy("Cypripedium calceolus",all=TRUE)
#' 
#' @export
wsl_taxonomy=function(sp_name = NULL, conf_match = 90, all = FALSE)
{
    # Search input name via GBIF backbone & error handling
    gbif.backbone = name_backbone(sp_name)

    if (gbif.backbone$matchType%in%"NONE") {
      cat("No species name found...","\n")
      return(NULL)

    } else if (gbif.backbone$confidence<conf_match) {
      cat("Confidence match not high enough...","\n")
      return(NULL)
    }

    # Extract key of accepted name
    if (gbif.backbone$status%in%"SYNONYM") {
     accep.key = gbif.backbone$acceptedUsageKey

    } else {
      accep.key = gbif.backbone$usageKey
    }

    # Extract accepted name and save it with its key in the prepared output
    accep.name = name_usage(accep.key,data="name")$data
    syn.syn = name_usage(accep.key,data="synonyms")$data

    # If missing codes, then we continue the search to find possible name correspondancy
    if (all) {

      # Extract children names
      syn.child = name_usage(accep.key,data="children")$data

      # Combine everything and search for related names (i.e. other string version)
      all.key = suppressWarnings(c(accep.name$key,syn.syn$key,syn.child$key))
      all.version = unlist(lapply(all.key,function(x){
        suppressWarnings(name_usage(x,data="related")$data$scientificName)
      }))

      # Extract all names
      a.n = suppressWarnings(accep.name$scientificName)
      s.n = suppressWarnings(syn.syn$scientificName)
      c.n = suppressWarnings(syn.child$scientificName)
      r.n = suppressWarnings(unique(all.version))
      all.names = c(a.n,s.n,c.n,r.n)

      # Extract all names and infos in a data.frame
      out = data.frame(Names = all.names,
        Status = c("ACCEPTED",
                 rep("SYNONYM",length(s.n)),
                 rep("CHILDREN",length(c.n)),
                 rep("RELATED",length(r.n))))
    
    } else {

      # Extract accepted names and synonyms
      out = data.frame(Names = suppressWarnings(c(accep.name$scientificName,syn.syn$scientificName)),
        Status = c("ACCEPTED",rep("SYNONYM",length(syn.syn$scientificName))))
    }
    return(out)
}

