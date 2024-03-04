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
get_status=function(sp_name = NULL, conf_match = 80, all = FALSE)
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
    iucn = try(rgbif::name_usage(accep.key,data="iucnRedListCategory")$data,silent=TRUE)
    if (methods::is(iucn,"try-error")) {iucn = "INTERNAL_GBIF_ERROR"} else {iucn = iucn$category}

     # Specific columns
    c.key = suppressWarnings(c(accep.key,syn.syn$key))
    c.sc = suppressWarnings(c(accep.name$scientificName,syn.syn$scientificName))
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
    out = data.frame(gbif_key = c.key,
      scientificName = c.sc,
      gbif_status = c.status,
      main.out,
      IUCN_status = iucn)

    return(out[!duplicated(out[,2]),])
}

