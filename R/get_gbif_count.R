### =========================================================================
### get_gbif_count
### =========================================================================
#' Count GBIF Occurrences Before Downloading Data
#'
#' Query GBIF and return the number of available records for a taxon using the
#' same name-matching logic as \code{get_gbif()}, but without downloading the
#' occurrence table.
#' 
#' @param sp_name Character string with the species name. Scientific names at
#' genus-species level are expected; fuzzy matching is available when
#' \code{search = FALSE}.
#' @param search Logical. If \code{TRUE} (default), use a strict GBIF backbone
#' search and keep only species-, subspecies-, or variety-level matches. If
#' \code{FALSE}, use a more permissive search and optionally rely on
#' \code{rank}, \code{phylum}, \code{class}, \code{order}, and
#' \code{family} to resolve ambiguous matches.
#' @param rank Character string giving the preferred rank to keep:
#' \code{"SPECIES"}, \code{"SUBSPECIES"}, or \code{"VARIETY"}. When
#' \code{NULL}, rank priority is inferred from \code{sp_name}.
#' @param phylum Optional phylum used to disambiguate alternative GBIF matches.
#' Particularly useful for hemihomonyms.
#' @param class Optional class used to disambiguate alternative GBIF matches.
#' @param order Optional order used to disambiguate alternative GBIF matches.
#' @param family Optional family used to disambiguate alternative GBIF matches.
#' @param conf_match Numeric confidence threshold between 0 and 100 for the
#' GBIF backbone match. Default is \code{80}.
#' @param geo Spatial object used to restrict the query extent. Accepted classes
#' are \code{Extent}, \code{SpatExtent}, \code{SpatialPolygon},
#' \code{SpatialPolygonDataFrame}, \code{SpatVector}, and \code{sf}. The
#' default \code{NULL} queries the whole globe.
#' @param has_xy Logical. If \code{TRUE} (default), count only records with
#' coordinates. If \code{FALSE}, count only records without coordinates. If
#' \code{NULL}, count all records.
#' @param spatial_issue Logical. If \code{FALSE} (default), count only records
#' without geospatial issues. If \code{TRUE}, count only records with
#' geospatial issues. If \code{NULL}, count all records.
#' @details The function mirrors the taxonomic matching strategy used by
#' \code{get_gbif()}, then reports both the total number of GBIF records and
#' the number retained after applying the chosen filters.
#' @return A numeric vector of length two giving the total number of GBIF
#' records and the number retained by the requested filters. A formatted summary
#' is also printed to the console.
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
					has_xy = TRUE,
          spatial_issue = FALSE) {
	
	######################################################
	### Stop messages
	######################################################


  	# General
  	check_character_vector(sp_name, "sp_name")
    check_logical(search, "search")
    check_numeric(conf_match, "conf_match")

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


  	# Set geo.ref for study area
    if (!is.null(geo) && !isTRUE(all.equal(geo, terra::ext()))) {
      if (class(geo)[1] %in% "sf") {
        geo <- terra::vect(geo)
      }
      if (!(class(geo) %in% "SpatExtent")) {
        geo <- terra::ext(geo)
      }
      geo.ref <- paste0(
        "POLYGON((", geo$xmin, " ", geo$ymin, ", ",
          geo$xmax, " ", geo$ymin, ", ",
          geo$xmax, " ", geo$ymax, ", ",
          geo$xmin, " ", geo$ymax, ", ",
          geo$xmin, " ", geo$ymin, "))"
      )
    } else {
      geo.ref <- NULL
    }


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

      # Filter by the supplied criteria if results are available
      if (!bsearch$matchType[1] %in% "NONE"){ 
        if (any(q.crit)){
          id.crit <- c("rank", "phylum", "class", "order", "family")[q.crit]
          p.crit <- unlist(list(rank, phylum, class, order, family)[q.crit])
          n.test <- id.crit %in% names(bsearch)

          if (any(n.test)){
            # selecting which
            id.crit2 <- id.crit[n.test]
            p.crit2 <- p.crit[n.test]

            # Apply the requested criteria
            for (i in seq_along(id.crit2)){
              bsearch <- bsearch[
                c(bsearch[, id.crit2[i]])[[1]] %in% p.crit2[i],
              ]

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
      
      # Continue with or without the supplied criteria
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

            # If only subspecies and variety are found, prioritize one
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

            if (any(bsearch$status %in% "ACCEPTED") &&
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

  	# Check number of records in total first
    gbif.total <- rgbif::occ_count(
      taxonKey = sp.key,
      hasCoordinate = NULL,
      hasGeospatialIssue = NULL,
      geometry = NULL
    )

    gbif.records <- rgbif::occ_count(
      taxonKey = sp.key,
      hasCoordinate = has_xy,
      hasGeospatialIssue = spatial_issue,
      geometry = geo.ref
    )

    # Print summary
    l <- sprintf(
      "%-29s : %10d",
      c("Total number (all records)", "Kept records"),
      c(gbif.total, gbif.records)
    )
    w <- max(nchar(l)) + 4
    cat("|",strrep("-",w-2),"|\n| ",
        paste(l,collapse=" |\n| ")," |\n|",strrep("-",w-2),"|\n",sep="")

    fmt <- function(x) if (is.null(x)) "NULL" else as.character(x)
    cat("| Kept records according to parameters:\n")

    # Print additional information depending if global or regional
    if (is.null(geo.ref)) {
      cat(sprintf("| spatial_issue = %s, has_xy = %s\n",
                  fmt(spatial_issue), fmt(has_xy)))
    } else {
      cat(sprintf(
        paste0(
          "| spatial_issue = %s, has_xy = TRUE by default ",
          "('geo' was set)\n"
        ),
        fmt(spatial_issue)
        )
      )
    }

  	# Cancel request if n==0
  	if (gbif.records == 0 ) {
  		cat("No species records found...","\n")
  		return(NULL)
  	}
  	return(c(gbif.total, gbif.records))
}
	
