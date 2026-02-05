### =========================================================================
### get_gbif_count
### =========================================================================
#' Retrieve GBIF observations count to select study case
#'
#' Follows the same code implementation as \code{get_gbif} for retrieving
#' the number of occurrences per species without downloading the data.
#' 
#' @param sp_name Character. Species name from which the user wants to retrieve
#' all existing GBIF names with associated taxonomy and IUCN status.
#' @param search Logical. If \code{TRUE} (default), the function will strictly
#' look for the most relevant result across the names given by \code{rgbif}
#' (only species, subspecies and variety allowed here), and give an error if
#' name matching is impeded by synonym duplicates. If \code{FALSE}, the
#' function will simply pick the first most relevant name (taxa rank higher
#' than species are allowed here). Also, unlike \code{search = TRUE}, fuzzy
#' search (~approximative name matching) is here allowed, and the \code{rank},
#' \code{phylum}, \code{class}, \code{order} and \code{family} parameters are
#' optionally used only if no convincing match is found. \code{FALSE} is
#' particularly useful if the given species name already includes the author.
#' @param rank Character. "SPECIES", "SUBSPECIES" or "VARIETY". If \code{NULL}
#' (default), the order of priority is (1) species, (2) subspecies and (3)
#' variety unless "subsp." or "var." is found in the \code{sp_name} parameter.
#' @param phylum Character (optional). What is the species' Phylum? Adds a
#' criteria to deal with alternative name matches and select the right synonym.
#' Available options are the GBIF Phylums
#' ((\href{https://www.gbif.org/species/search}{listed per Kingdom/Phylum})).
#' If \code{search = FALSE}, only used if no direct match is found.
#' If \code{search = FALSE}, only used if no direct match is found.
#' @param class Character (optional). What is the species' Class? Same as above
#' but at the finer class level. Available options are the GBIF Classes
#' (same url). If \code{search = FALSE}, only used if no direct match is found.
#' @param order Character (optional). What is the species' Order? Same as above
#' but at the finer order level. Available options are the GBIF Orders (same
#' url). If \code{search = FALSE}, only used if no direct match is found.
#' @param family Character (optional). What is the species' Family? Same as
#' above but at the finer family level. Available options are the GBIF Families
#' (same url). If \code{search = FALSE}, used only if no direct match is found.
#' @param conf_match Numeric from 0 to 100. Determine the confidence threshold
#' of match between \code{sp_name} and the GBIF backbone taxonomy.
#' Default is 90.
#' @param geo Object of class \code{Extent}, \code{SpatExtent},
#' \code{SpatialPolygon}, \code{SpatialPolygonDataframe}, \code{SpatVector}
#' or \code{sf} (WGS84) to define the study's area extent. Default is
#' \code{NULL}, i.e., the whole globe.
#' @param no_xy Logical. Only number of records with coordinates are given.
#' Default is \code{FALSE}. If \code{TRUE}, number of records without
#' coordinates is given.
#' @details Implements the same search result when
#' (\href{https://www.gbif.org}{GBIF}) is employed, i.e., based on the
#' input taxa name, all species records related to its accepted name
#' and synonyms are extracted.
#' @return A simple print of the number of observations found.
#' @references
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the
#' global biodiversity information facility API. 10.5281/zenodo.6023735
#' @example inst/examples/get_gbif_count_help.R
#' @importFrom terra ext vect
#' @importFrom rgbif name_backbone occ_search
#' @importFrom CoordinateCleaner cd_ddmm cd_round
#' @export
get_gbif_count <- function(sp_name = NULL,
					search = TRUE,
					rank = NULL,
					phylum = NULL,
					class = NULL,
					order = NULL,
					family = NULL,
					conf_match = 80,
					geo = NULL,
					no_xy = FALSE) {
	
	######################################################
	### Stop messages
	######################################################


	# General
	check_character_vector(sp_name, "sp_name")
  check_logical(search, "search")
  check_numeric(conf_match, "conf_match")
  check_logical(no_xy, "no_xy")

  # Geo class
  spatial.class <- c("Extent", "SpatExtent", "SpatialPolygon",
  	"SpatialPolygonDataframe", "SpatVector", "sf")
  if (!is.null(geo)) {
  	 if (!class(geo)[1] %in% spatial.class) {
     	stop("Wrong 'geo' class (not a spatial object)...")
     }
  }

	
	#############################################################
	### Convert to SpatVect if class is 'sf' + convert to extent
	#############################################################


	# For study area
	if (!is.null(geo)) {
		if (class(geo)[1] %in% "sf") {geo <- terra::vect(geo)}
		if (!(class(geo) %in% "SpatExtent")) {geo <- terra::ext(geo)}

	} else {
		geo <- terra::ext(-180, 180, -90, 90)
	}


	######################################################
	### Parameter config
	######################################################


	# Set geo.ref
	geo.ref <- paste0(
		"POLYGON((", geo$xmin, " ", geo$ymin, ", ",
			geo$xmax, " ", geo$ymin, ", ",
			geo$xmax, " ", geo$ymax, ", ",
			geo$xmin, " ", geo$ymax, ", ",
			geo$xmin, " ", geo$ymin, "))"
	)


	######################################################
	### Backbone harmonization and N records
	######################################################


	if (!search){
    # Search input name via fuzzy match and direct search
    bsearch <- rgbif::name_backbone(sp_name,
    	rank = rank,
    	phylum = phylum,
    	class = class,
    	order = order,
    	family = family,
    	verbose = FALSE,
    	strict = FALSE)

  } else {
    # Search input name via strict match and refined search
    bsearch <- rgbif::name_backbone(sp_name,
    	verbose = TRUE,
    	strict = TRUE)

    q.crit <- !vapply(
    	list(rank, phylum, class, order, family),
    	is.null,
    	logical(1)
    )

    # Filter by given criterias if results
    if (!bsearch$matchType[1] %in% "NONE"){ 
      if (any(q.crit)){
        id.crit <- c("rank", "phylum", "class", "order", "family")[q.crit]
        p.crit <- unlist(list(rank, phylum, class, order, family)[q.crit])
        n.test <- id.crit %in% names(bsearch)

        if (any(n.test)){
          # selecting which
          id.crit2 <- id.crit[n.test]
          p.crit2 <- p.crit[n.test]

          # Apply the rigth criterias
          for (i in seq_along(id.crit2)){
            bsearch <- bsearch[c(bsearch[, id.crit2[i]])[[1]] %in% p.crit2[i], ]

            if (nrow(bsearch) == 0){
              bsearch <- data.frame(matchType = "NONE")
            }
          }
        }

        if (!all(n.test)){
          pp <- paste(id.crit[!n.test], collapse = ", ")
          warning(
          		"'", pp, 
          		"' level(s) not available in GBIF, could not be employed...",
          		"\n"
          )
        }
      }
    }
    
    # Normal procedure with or without criterias
    if (nrow(bsearch) > 1){
      if (all(!bsearch$rank %in% c("SPECIES", "SUBSPECIES", "VARIETY"))){
        cat("Not match found...", "\n")
        return(NULL)

      } else {
        s.keep <- bsearch[bsearch$rank %in%
        								c("SPECIES" ,"SUBSPECIES" ,"VARIETY"),]
        s.keep <- s.keep[s.keep$status %in% c("ACCEPTED", "SYNONYM"),]
        if (nrow(s.keep) == 0){
          cat("Not match found...", "\n")
          return(NULL)

        } else if (nrow(s.keep) > 1){

          # If we only find subpsecies and variety, we need to prioritize
          if (all(s.keep$rank %in% c("VARIETY", "SUBSPECIES"))){
            if ("var." %in% strsplit(sp_name," ")[[1]]){
              bsearch <- s.keep[s.keep$rank %in% "VARIETY", ]
              if (nrow(bsearch) == 0){
                bsearch <- s.keep[s.keep$rank %in% "SUBSPECIES", ]
              }

            } else {
              bsearch <- s.keep[s.keep$rank %in% "SUBSPECIES", ]
            }

          } else {
            bsearch <- s.keep[s.keep$rank %in% "SPECIES", ]
          }

          focp.key <- c("familyKey", "orderKey", "classKey", "phylumKey")
          coltax <- focp.key  %in% colnames(bsearch)
          key.test <- bsearch[ ,focp.key[coltax]]

          if (any(bsearch$status %in% "ACCEPTED") &
          						length(unique(key.test[, 1])) == 1){
          	bsearch <- bsearch[bsearch$status %in% "ACCEPTED", ]
          }

        } else {
          bsearch <- s.keep
        }
        # If not the same species overall return empty
        s.usp <- length(unique(bsearch$speciesKey)) == 1
        if (!s.usp){
          cat("No synonyms distinction could be made.",
          	"Consider using phylum/class/order/family...","\n")
          return(NULL)

        } else {
          bsearch <- bsearch[1, ]
        } 
      }
    }
  }

  if (bsearch$matchType %in% "NONE") {
    cat("No species name found...","\n")
    return(NULL)
  }

  if (bsearch$confidence[1] < conf_match) {
    cat("Confidence match not high enough...","\n")
    return(NULL)
  }  

	# Get the accepetedKey
	if (bsearch$status %in% "SYNONYM"){
		sp.key <- bsearch$acceptedUsageKey

	} else {
		sp.key <- bsearch$usageKey
	}

	# Check number of records in 'geo' first
	gbif.records <- rgbif::occ_count(taxonKey = sp.key,
									 hasCoordinate = !no_xy,
									 hasGeospatialIssue = FALSE,
									 geometry = geo.ref)

	cat("Total number of records:", gbif.records,"\n")

	# Cancel request if n==0
	if (gbif.records == 0 ) {
		cat("No species records found...","\n")
		return(NULL)
	}
	return(gbif.records)
}
	