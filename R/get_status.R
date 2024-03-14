### ==================================================================
### get_status
### ==================================================================
#' Retrieve from GBIF the IUCN and taxonomy status of a specific Taxa
#'
#' Generates, based on a given species name, its IUCN red list status and a list of
#' all scientific names (accepted, synonyms) found in the GBIF backbone taxonomy and used to
#' download the data in get_gbif(). Children and related doubtful names not used to download
#' the data may also be extracted.
#'
#' @param sp_name Character. Species name from which the user wants to retrieve all existing GBIF names.
#' @param rank Character. "SPECIES", "SUBSPECIES" or "VARIETY". If NULL (default), the order of priority
#' is (1) species, (2) subspecies and (3) variety unless "subsp." or "var." is found in 'sp_name'.
#' @param phylum Character. Optional. What is the species' Phylum? Adds a criteria to deal with alternative
#' name matches and select the right synonym. Available options are the GBIF Phylums
#' (listed per Kingdom --> https://www.gbif.org/species/1).
#' @param class Character. Optional. What is the species' Class? Same as above but at the finer class level.
#' Available options are the GBIF Classes (same url).
#' @param order Character. Optional. What is the species' Order? Same as above but at the finer order level.
#' Available options are the GBIF Orders (same url).
#' @param family Character. Optional. What is the species' Family? Same as above but at the finer order level.
#' Available options are the GBIF Orders (same url).
#' @param conf_match Numeric. From 0 to 100. Determine the confidence
#' threshold of match of 'sp_name' with the GBIF backbone taxonomy. Default is 90.
#' @param all Logical. Default is FALSE. Should all species names be retrieved or only
#' the accepted name and its synonyms?
#' 
#' @return Data.frame with three columns: (1) GBIF taxonomic key, (2) scientificName and
#' (3) Backbone Taxonomy Status, (4) Genus, (5) Family, (6) Order, (7) Phylum and (8) IUCN status
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches on how to retrieve
#' scientific names from the GBIF backbone taxonomy.
#' @examples
#' get_status("Cypripedium calceolus",all=FALSE)
#' get_status("Cypripedium calceolus",all=TRUE)
#' 
#' @export
#' @importFrom rgbif name_backbone name_usage
#' @importFrom methods is
get_status=function(sp_name = NULL,
                    rank = NULL,
                    phylum = NULL,
                    class = NULL,
                    order = NULL,
                    family = NULL,
                    conf_match = 80,
                    all = FALSE)
{
    # Search input name via GBIF backbone & error handling
    bone.search = rgbif::name_backbone(sp_name,
                                       rank = rank,
                                       phylum = phylum,
                                       class = class,
                                       order = order,
                                       family = family,
                                       verbose = TRUE,
                                       strict = TRUE)

    if (!is.null(c(rank,phylum,class,order,family))){
      bone.search = bone.search[1,]
    
    } else {
      if (nrow(bone.search)>1){
        if (all(!bone.search$rank%in%c("SPECIES","SUBSPECIES","VARIETY"))){
          cat("Not match found...","\n")
          return(NULL)

        } else {
          s.keep = bone.search[bone.search$rank%in%c("SPECIES","SUBSPECIES","VARIETY"),]
          s.keep = s.keep[s.keep$status%in%c("ACCEPTED","SYNONYM"),]
          s.keep = s.keep[s.keep$matchType%in%"EXACT",]
          if (nrow(s.keep)==0){
            cat("Not match found...","\n")
            return(NULL)

          } else if (nrow(s.keep)>1){

            # If we only find subpsecies and variety, we need to (default) prioritize
            if (all(s.keep$rank%in%c("VARIETY","SUBSPECIES"))){
              if ("var."%in%strsplit(sp_name," ")[[1]]){
                bone.search = s.keep[s.keep$rank%in%"VARIETY",]
                if (nrow(bone.search)==0){
                  bone.search = s.keep[s.keep$rank%in%"SUBSPECIES",]
                }

              } else {
                bone.search = s.keep[s.keep$rank%in%"SUBSPECIES",]
              }

            } else {
              bone.search = s.keep[s.keep$rank%in%"SPECIES",]
            }

            if (any(bone.search$status%in%"ACCEPTED") & length(unique(bone.search$family))==1){
              bone.search = bone.search[bone.search$status%in%"ACCEPTED",]
            }

          } else {
            bone.search = s.keep
          }
          # If not the same species overall return NULL
          s.usp = length(unique(bone.search$species))==1
          if (!s.usp){
            cat("No synonyms distinction could be made. Consider using phylum/class/order/family...","\n")
            return(NULL)

          } else {
            bone.search = bone.search[1,]
          } 
        }
      }
    }

    if (bone.search$matchType%in%"NONE") {
      cat("No species name found...","\n")
      return(NULL)
    }

    if (bone.search$confidence[1]<conf_match) {
      cat("Confidence match not high enough...","\n")
      return(NULL)
    }  

    # Extract key of accepted name
    if (bone.search$status%in%"SYNONYM") {
      accep.key = bone.search$acceptedUsageKey
    } else {
      accep.key = bone.search$usageKey
    }

    # Extract accepted name and save it with its key in the prepared output
    accep.name = rgbif::name_usage(accep.key,data="name")$data
    syn.syn = rgbif::name_usage(accep.key,data="synonyms")$data
    main.dat =  rgbif::name_usage(accep.key,data="all")$data
    iucn = try(rgbif::name_usage(accep.key,data="iucnRedListCategory")$data,silent=TRUE)
    if (methods::is(iucn,"try-error")) {iucn = "NOT_FOUND"} else {iucn = iucn$category}

     # Specific columns
    c.key = suppressWarnings(c(accep.key,syn.syn$key))
    c.sc = suppressWarnings(c(accep.name$scientificName,syn.syn$scientificName))
    c.can = suppressWarnings(c(accep.name$canonicalName,syn.syn$canonicalName))
    c.status = c("ACCEPTED",rep("SYNONYM",length(suppressWarnings(syn.syn$scientificName))))

    # If missing codes, then we continue the search to find possible name correspondence
    if (all) {

      # Combine everything and search for related names (i.e. other string version)
      all.key = suppressWarnings(c(accep.key,syn.syn$key,main.dat$key))
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
      c.n = suppressWarnings(main.dat[,c("key","scientificName")])
      r.n = suppressWarnings(unique(do.call("rbind",all.version)))

      # Conditions for synonymy
      s.n = try(suppressWarnings(syn.syn[,c("key","scientificName")]),silent=TRUE)
      if (class(s.n)[1]%in%"try-error") {s.n = data.frame(key=NULL,scientificName=NULL)}
      all.names = rbind(a.n,s.n,c.n,r.n)

      # Specific columns
      c.key = all.names$key
      c.sc = suppressWarnings(all.names$scientificName)
      c.can = suppressWarnings(all.names$canonicalName)
      c.status = c("ACCEPTED",
                rep("SYNONYM",nrow(s.n)),
                rep("CHILDREN",nrow(c.n)),
                rep("RELATED",nrow(r.n)))
    }

    # Which is null in main.dat for Genus, Family, Order, Phyllum?
    exist.not = c("genus","order","family","phylum")%in%names(main.dat)
    main.out = data.frame(Genus="Unknown",Family="Unknown",Order="Unknown",Phylum="Unknown")
    main.out[exist.not] = main.dat[,c("genus","order","family","phylum")[exist.not]]

    # Extract accepted names and synonyms
    out = data.frame(canonicalName = c.can,
      gbif_key = c.key,
      scientificName = c.sc,
      gbif_status = c.status,
      main.out,
      IUCN_status = iucn)

    return(out[!duplicated(out[,2]),])
}

