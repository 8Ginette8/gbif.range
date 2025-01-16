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
#' @param sp_name Character. Species name from which the user wants to retrieve all existing GBIF names
#' with associated taxonomy and IUCN status.
#' @param search Logical. If TRUE, the function will strictly look for the most relevant result, based
#' on a list of names given by rgbif (only species, subspecies and variety allowed here), and give an error
#' if name matching was impeded by synonym duplicates. If FALSE, the function will simply pick the first most
#' relevant name from the list (higher taxa level than species allowed here). Also, unlike search=TRUE,
#' fuzzy search (~approximative name match) is here allowed, and the 'rank', phylum', 'class', order'
#' and 'family' parameters are optionally used only if no convincing name match is found. FALSE is
#' particularly useful if the given species name already include the author.
#' @param rank Character. "SPECIES", "SUBSPECIES" or "VARIETY". If NULL (default), the order of priority
#' is (1) species, (2) subspecies and (3) variety unless "subsp." or "var." is found in 'sp_name'.
#' @param phylum Character. Optional. What is the species' Phylum? Adds a criteria to deal with alternative
#' name matches and select the right synonym. Available options are the GBIF Phylums
#' (listed per Kingdom --> https://www.gbif.org/species/1). If search = FALSE, used only if no direct match
#' is found.
#' @param class Character. Optional. What is the species' Class? Same as above but at the finer class level.
#' Available options are the GBIF Classes (same url). If search = FALSE, used only if no direct match
#' is found.
#' @param order Character. Optional. What is the species' Order? Same as above but at the finer order level.
#' Available options are the GBIF Orders (same url). If search = FALSE, used only if no direct match
#' is found.
#' @param family Character. Optional. What is the species' Family? Same as above but at the finer family level.
#' Available options are the GBIF Orders (same url). If search = FALSE, used only if no direct match
#' is found.
#' @param conf_match Numeric. From 0 to 100. Determine the confidence threshold of match of 'sp_name' with
#' the GBIF backbone taxonomy. Default is 90.
#' @param all Logical. Default is FALSE. Should all species names be retrieved or only the accepted name
#' and its synonyms?
#' 
#' @return Data.frame with nine columns: (0) Simplified name, (1) GBIF taxonomic key, (2) scientificName, 
#' (3) Backbone Taxonomy Status, (4) Genus, (5) Family, (6) Order, (7) Phylum, (8) IUCN status and
#' sp_nameMatch informing how well the input sp_name has matched with the found synonyms.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' @seealso The rgbif package for additional and more general approaches on how to retrieve
#' scientific names from the GBIF backbone taxonomy.
#' @example inst/examples/get_status_help.R
#' @importFrom rgbif name_backbone name_usage
#' @importFrom methods is
#' @export
get_status=function(sp_name = NULL,
                    search = TRUE,
                    rank = NULL,
                    phylum = NULL,
                    class = NULL,
                    order = NULL,
                    family = NULL,
                    conf_match = 80,
                    all = FALSE)
{
  # Error message
  if (is.na(sp_name)|is.null(sp_name)){
    stop("Given 'sp_name' is NA or NULL, a string must be provided...")
  }

  # Empty output
  e_output = data.frame(canonicalName = NA,
                       rank = NA,
                       gbif_key = NA,
                       scientificName = NA,
                       gbif_status = NA,
                       Genus = NA,
                       Family = NA,
                       Order = NA,
                       Class = NA,
                       Phylum = NA,
                       IUCN_status = NA,
                       sp_nameMatch = NA)

  # Search
  if (!search){
    # Search input name via fuzzy match and direct search
    bone_search = rgbif::name_backbone(sp_name,
                                     rank = rank,
                                     phylum = phylum,
                                     class = class,
                                     order = order,
                                     family = family,
                                     verbose = FALSE,
                                     strict = FALSE)
  } else {
    # Search input name via strict match and refined search
    bone_search = rgbif::name_backbone(sp_name,
                                       verbose = TRUE,
                                       strict = TRUE)

    q_crit = !sapply(list(rank,phylum,class,order,family),is.null) 

    # Filter by given criterias if results
    if (!bone_search$matchType[1]%in%"NONE"){ 
      if (any(q_crit)){
        id_crit = c("rank","phylum","class","order","family")[q_crit]
        p_crit = unlist(list(rank,phylum,class,order,family)[q_crit])
        n_test = id_crit%in%names(bone_search)
        if (any(n_test)){
          # selecting which
          id_crit2 = id_crit[n_test]
          p_crit2 = p_crit[n_test]
          # Apply the rigth criterias
          for (i in 1:length(id_crit2)){
            bone_search = bone_search[c(bone_search[,id_crit2[i]])[[1]]%in%p_crit2[i],]
            if (nrow(bone_search)==0){
              bone_search = data.frame(matchType="NONE")
            }
          }
        }
        if (!all(n_test)){
          pp = paste(id_crit[!n_test],collapse=", ")
          warning(paste0("'",pp,"' level(s) not available for this taxa in GBIF, could not be employed..."))
        }
      }
    }
    
    # Normal procedure with or without criterias
    if (nrow(bone_search)>1){
      if (all(!bone_search$rank%in%c("SPECIES","SUBSPECIES","VARIETY"))){
        cat("Not match found...","\n")
        return(e_output)

      } else {
        s_keep = bone_search[bone_search$rank%in%c("SPECIES","SUBSPECIES","VARIETY"),]
        s_keep = s_keep[s_keep$status%in%c("ACCEPTED","SYNONYM"),]
        if (nrow(s_keep)==0){
          cat("Not match found...","\n")
          return(e_output)

        } else if (nrow(s_keep)>1){

          # If we only find subpsecies and variety, we need to (default) prioritize
          if (all(s_keep$rank%in%c("VARIETY","SUBSPECIES"))){
            if ("var."%in%strsplit(sp_name," ")[[1]]){
              bone_search = s_keep[s_keep$rank%in%"VARIETY",]
              if (nrow(bone_search)==0){
                bone_search = s_keep[s_keep$rank%in%"SUBSPECIES",]
              }

            } else {
              bone_search = s_keep[s_keep$rank%in%"SUBSPECIES",]
            }

          } else {
            bone_search = s_keep[s_keep$rank%in%"SPECIES",]
          }

          coltax = c("familyKey","orderKey","classKey","phylumKey")%in%colnames(bone_search)
          key_test = bone_search[,c("familyKey","orderKey","classKey","phylumKey")[coltax]]

          if (any(bone_search$status%in%"ACCEPTED") & length(unique(key_test[,1][[1]]))==1){
            bone_search = bone_search[bone_search$status%in%"ACCEPTED",]
          }

        } else {
          bone_search = s_keep
        }
        # If not the same species overall return NULL
        s_usp = length(unique(bone_search$speciesKey))==1
        if (!s_usp){
          cat("No synonyms distinction could be made. Consider using phylum/class/order/family...","\n")
          return(e_output)

        } else {
          bone_search = bone_search[1,]
        } 
      }
    }
  }

  if (bone_search$matchType%in%"NONE") {
    cat("No species name found...","\n")
    return(e_output)
  }

  if (bone_search$confidence[1]<conf_match) {
    cat("Confidence match not high enough...","\n")
    return(e_output)
  }  

   # Keep original scientific name
  sc_name = suppressWarnings(bone_search$scientificName)

  # Extract key of accepted name
  if (bone_search$status%in%"SYNONYM") {
    accep_key = bone_search$acceptedUsageKey
  } else {
    accep_key = bone_search$usageKey
  }

  # Extract accepted name and save it with its key in the prepared output
  accep_name = rgbif::name_usage(accep_key,data="name")$data
  syn_syn = rgbif::name_usage(accep_key,data="synonyms")$data
  main_dat =  rgbif::name_usage(accep_key,data="all")$data
  iucn = try(rgbif::name_usage(accep_key,data="iucnRedListCategory")$data,silent=TRUE)
  if (methods::is(iucn,"try-error")) {iucn = "NOT_FOUND"} else {iucn = iucn$category}

  # Avoid the canonicalName error (sometimes the column is absent wtf...)
  if (!"canonicalName"%in%colnames(main_dat)){
    accep_name$canonicalName = NA
    syn_syn$canonicalName = NA
    main_dat$canonicalName = NA
  }

  # Specific columns
  c_key = suppressWarnings(c(accep_key,syn_syn$key))
  c_sc = suppressWarnings(c(accep_name$scientificName,syn_syn$scientificName))
  c_can = suppressWarnings(c(accep_name$canonicalName,syn_syn$canonicalName))
  c_status = c("ACCEPTED",rep("SYNONYM",length(suppressWarnings(syn_syn$scientificName))))

  # If all=TRUE, then we continue the search to find possible name correspondence
  if (all) {

    # Combine everything and search for related names (i.e. other string version)
    all_key = suppressWarnings(c(accep_key,syn_syn$key,main_dat$key))
    all_version = lapply(all_key,function(x){
      out = suppressWarnings(rgbif::name_usage(x,data="related")$data)
      if (is.null(out)|nrow(out)==0){
        return(NULL)
      } else {
        n_ref = c("canonicalName","scientificName")
        r_col = n_ref[n_ref%in%names(out)]
        r_out = data.frame(canonicalName=rep(NA,nrow(out)),scientificName=rep(NA,nrow(out)))
        r_out[,r_col] = out[,r_col]
        return(data.frame(key = x,
                          canonicalName = r_out$canonicalName,
                          scientificName = r_out$scientificName))
      }
    })

    # Extract all names
    accep_n = suppressWarnings(accep_name[,c("canonicalName","key","scientificName")])
    accep_n$key = accep_key
    c_n = suppressWarnings(main_dat[,c("canonicalName","key","scientificName")])
    r_n = suppressWarnings(unique(do.call("rbind",all_version)))

    # Conditions for synonymy
    syn_n = try(suppressWarnings(syn_syn[,c("canonicalName","key","scientificName")]),silent=TRUE)
    if (class(syn_n)[1]%in%"try-error") {syn_n = data.frame(canonicalName=NULL,key=NULL,scientificName=NULL)}
    all_names = rbind(accep_n,syn_n,c_n,r_n)

    # Specific columns
    c_key = all_names$key
    c_sc = suppressWarnings(all_names$scientificName)
    c_can = suppressWarnings(all_names$canonicalName)
    c_status = c("ACCEPTED",
              rep("SYNONYM",nrow(syn_n)),
              rep("CHILDREN",nrow(c_n)),
              rep("RELATED",nrow(r_n)))
  }

  # Which is null in main.dat for Genus, Family, Order, Phyllum?
  exist_not = c("genus","family","order","class","phylum")%in%names(main_dat)
  main_out = data.frame(Genus=NA,Family=NA,Order=NA,Class=NA,Phylum=NA)
  main_out[exist_not] = main_dat[,c("genus","family","order","class","phylum")[exist_not]]

  # Extract accepted names and synonyms
  e_output = data.frame(canonicalName = c_can,
                        rank = bone_search$rank,
                        gbif_key = c_key,
                        scientificName = c_sc,
                        gbif_status = c_status,
                        main_out,
                        IUCN_status = iucn,
                        sp_nameMatch = bone_search$matchType)

  # Remove duplicated and change sp_nameMatch to know which row related to the sp_name input
  e_output = e_output[!duplicated(e_output$scientificName),]
  e_output[e_output$scientificName%in%sc_name,"sp_nameMatch"][1] = "INPUT"

  return(e_output)
}
