### ==================================================================
### get_taxonomy
### ==================================================================
#' Retrieve from GBIF all scientific names of a specific Taxa
#'
#' Generates, based on a given species name, a list of all its scientific names
#' (accepted, synonyms) found in the GBIF backbone taxonomy and used to download the data in
#' get_gbif(). Children and related doubtful names not used to download the data may also be extracted.
#'
#' @param sp_name Character. Species name from which the user wants to retrieve all existing GBIF names.
#' @param conf_match Numeric. From 0 to 100. Determine the confidence
#' threshold of match of 'sp_name' with the GBIF backbone taxonomy. Default is 90.
#' @param all Logical. Default is FALSE. Should all species names be retrieved or only
#' the accepted name and its synonyms?
#' 
#' @return Data.frame with three columns: (1) GBIF taxonomic key, (2) scientificName and
#' (3) Backbone Taxonomy Status.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches on how to retrieve
#' scientific names from the GBIF backbone taxonomy.
#' @examples
#' get_taxonomy("Cypripedium calceolus",all=FALSE)
#' get_taxonomy("Cypripedium calceolus",all=TRUE)
#' 
#' @export
get_taxonomy=function(sp_name = NULL, conf_match = 80, all = FALSE)
{
    # Search input name via GBIF backbone & error handling
    gbif.backbone = rgbif::name_backbone(sp_name)

    if (gbif.backbone$matchType%in%"NONE") {
      cat("No species name found...","\n")
      return(NULL)

    } else if (gbif.backbone$confidence<conf_match) {
      cat("Confidence match not high enough...","\n")
      return(NULL)
    
    } else if (!gbif.backbone$rank%in%c("SPECIES","SUBSPECIES","VARIETY")) {
      gbif.backbone = rgbif::name_backbone(sp_name,verbose=TRUE)[2,]
      if (gbif.backbone$confidence < conf_match | is.na(gbif.backbone$status)) {
        cat("Not match found...","\n")
        return(NULL)
      }
    }

    # Extract key of accepted name
    if (gbif.backbone$status%in%"SYNONYM") {
     accep.key = gbif.backbone$acceptedUsageKey

    } else {
      accep.key = gbif.backbone$usageKey
    }

    # Extract accepted name and save it with its key in the prepared output
    accep.name = rgbif::name_usage(accep.key,data="name")$data
    syn.syn = rgbif::name_usage(accep.key,data="synonyms")$data
    main.dat =  rgbif::name_usage(accep.key,data="all")$data

    # If missing codes, then we continue the search to find possible name correspondence
    if (all) {

       # Extract children names
      syn.child =  main.dat$data

      # Combine everything and search for related names (i.e. other string version)
      all.key = suppressWarnings(c(accep.key,syn.syn$key,syn.child$key))
      all.version = lapply(all.key,function(x){
        out = suppressWarnings(rgbif::name_usage(x,data="related")$data$scientificName)
        if (is.null(out)){
          return(NULL)
        } else {
          return(data.frame(key=x,scientificName=out))
        }
      })

      # Extract all names
      a.n = suppressWarnings(accep.name[,c("key","scientificName")])
      a.n$key = accep.key
      s.n = suppressWarnings(syn.syn[,c("key","scientificName")])
      c.n = suppressWarnings(syn.child[,c("key","scientificName")])
      r.n = suppressWarnings(unique(do.call("rbind",all.version)))
      all.names = rbind(a.n,s.n,c.n,r.n)

      # Extract all names and infos in a data.frame
      out = data.frame(key=all.names$key,
        scientificName=suppressWarnings(all.names$scientificName),
        status=c("ACCEPTED",
                rep("SYNONYM",nrow(s.n)),
                rep("CHILDREN",nrow(c.n)),
                rep("RELATED",nrow(r.n))))

      return(out[!duplicated(out[,2]),])
    
    } else {

      # Extract accepted names and synonyms
      out = data.frame(key=suppressWarnings(c(accep.key,syn.syn$key)),
        scientificName=suppressWarnings(c(accep.name$scientificName,syn.syn$scientificName)),
        status=c("ACCEPTED",rep("SYNONYM",length(suppressWarnings(syn.syn$scientificName)))),
        genus=main.dat$genus,
        family=main.dat$family,
        order=main.dat$order,
        phylum=main.dat$phylum)
    }
    
    return(out[!duplicated(out[,2]),])
}

