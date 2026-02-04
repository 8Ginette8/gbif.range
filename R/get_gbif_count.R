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
#' (listed per Kingdom/Phylum: https://www.gbif.org/species/search).
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
#' @param no_xy Logical. Only records with coordinates are downloaded.
#' Default is \code{FALSE}. If \code{TRUE}, records with no coordinates are
#' also downloaded.
#' @param add_infos Character. Infos that may be added to the default output
#' information. Default IDs contain "taxonKey",
#' "scientificName", "acceptedTaxonKey", "acceptedScientificName",
#' "individualCount", "decimalLatitude", "decimalLongitude", "basisOfRecord",
#' "coordinateUncertaintyInMeters", "countryCode", "country", "year",
#' "datasetKey", "institutionCode", "publishingOrgKey", "taxonomicStatus",
#' "taxonRank" and "degreeOfEstablishment".
#' List of IDs may be found at:
#' https://www.gbif.org/developer/occurrence.
#' @param time_period Numerical vector. Observations will be downloaded
#' according to the chosen year range. Default is \code{c(1000,3000)}.
#' Observations with \code{NA} are kept by default.
#' @param identic_xy Logical. Should records with identical xy be kept?
#' Default is \code{FALSE}.
#' @param wConverted_xy Logical. Should incorrectly converted lon/lat be
#' kept? Default is \code{TRUE}. Otherwise, implements an approximate version
#' of \code{cd_ddmm()} from the \code{CoordinateCleaner} R package. See this
#' package for more advanced options.
#' @param centroids Logical. Should species records from raster centroids be
#' kept? Default is \code{TRUE}. Uses \code{cd_round()} from the
#' \code{CoordinateCleaner} R package.
#' @param ntries Numeric. In case of internal errors (GBIF server or
#' \code{rgbif} R package), how many download attempts should
#' \code{get_gbif()} request? Default is \code{10} with a 2 seconds interval
#' between tries. If the attempts failed, an empty data.frame is return by
#' default.
#' @param error_skip Logical. Should the search process continues if
#' \code{ntries} failed ?
#' @param occ_samp Numeric. Determine how many GBIF occurrences will be
#' sampled per geographic tiles of the fragmented study area. Default is the
#' maximum number of GBIF observations found in a tile (i.e. ~10'000 records).
#' A lower number may be set (< 10'000) if the user only wants a sample of the
#' species GBIF observations, hence increasing the download process.
#' @param should_use_occ_download Logical. If \code{TRUE}, \code{get_gbif()}
#' will use the \code{rgbif::occ_download()} instead of
#´\code{rgbif::occ_search()}.  This requires GBIF credentials! Defaults to
#' \code{FALSE}.
#' @param occ_download_user Character. GBIF username.
#' Required if \code{should_use_occ_download = TRUE}.
#' @param occ_download_pwd Character. GBIF password.
#' Required if \code{should_use_occ_download = TRUE}.
#' @param occ_download_email Character. GBIF email.
#' Required if \code{should_use_occ_download = TRUE}.
#' @param ... Additonnal parameters for the function \code{cd_round()} of
#' the \code{CoordinateCleaner} R package.
#' @details Implements the same search result when
#' www.gbif.org is employed, i.e., based on the input taxa name, all species
#' records related to its accepted name and synonyms are extracted.
#' @return A simple print of the number of observations found.
#' @references
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the
#' global biodiversity information facility API. 10.5281/zenodo.6023735


#' @example inst/examples/get_gbif_count_help.R
#' @importFrom terra ext vect
#' @importFrom rgbif name_backbone occ_search
#' @importFrom CoordinateCleaner cd_ddmm cd_round
#' @export
get_gbif <- function(sp_name = NULL,
					search = TRUE,
					rank = NULL,
					phylum = NULL,
					class = NULL,
					order = NULL,
					family = NULL,
					conf_match = 80,
					geo = NULL,
					no_xy = FALSE,
					add_infos = NULL,
					time_period = c(1000,3000),
					identic_xy = FALSE,
					wConverted_xy = TRUE,
					centroids = FALSE,
					ntries = 10,
					error_skip = TRUE,
					occ_samp = 10000,
					should_use_occ_download = FALSE,
					occ_download_user = NULL,
					occ_download_pwd = NULL,
					occ_download_email = NULL,
					...) {
	
	######################################################
	### Stop messages
	######################################################


	# General
	check_character_vector(sp_name, "sp_name")
  check_logical(search, "search")
  check_numeric(conf_match, "conf_match")
  check_logical(no_xy, "no_xy")
  check_numeric_range(time_period, "time_period", 2)
  check_logical(identic_xy, "identic_xy")
  check_logical(wConverted_xy, "wConverted_xy")
  check_logical(centroids, "centroids")
  check_numeric(ntries, "ntries")
  check_logical(error_skip, "error_skip")
  check_numeric(occ_samp, "occ_samp")
  check_logical(should_use_occ_download, "should_use_occ_download")

  # Add_infos
  if (all(!is.null(add_infos))) {
  	if (all(is.na(add_infos)))
  	stop("NA values not allowed for given 'add_infos'")
  }

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
        return(e.output)

      } else {
        s.keep <- bsearch[bsearch$rank %in%
        								c("SPECIES" ,"SUBSPECIES" ,"VARIETY"),]
        s.keep <- s.keep[s.keep$status %in% c("ACCEPTED", "SYNONYM"),]
        if (nrow(s.keep) == 0){
          cat("Not match found...", "\n")
          return(e.output)

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
          return(e.output)

        } else {
          bsearch <- bsearch[1, ]
        } 
      }
    }
  }

  if (bsearch$matchType %in% "NONE") {
    cat("No species name found...","\n")
    return(e.output)
  }

  if (bsearch$confidence[1] < conf_match) {
    cat("Confidence match not high enough...","\n")
    return(e.output)
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

	cat(">>>>>>>> Total number of records:", gbif.records,"\n")

	# Cancel request if n==0
	if (gbif.records == 0 ) {
		cat("No species records found...","\n")
		return(e.output)
	}


	######################################################
	### Find the tiles where obs. < 10'000 observations
	######################################################


	## 1) If species records > 10'000, search for the optimum tiles
	if (!should_use_occ_download && gbif.records > 10000) {
		# Start with 10 tiles
		tile.100 <- make_tiles(geo, ntiles = 10, sext = TRUE)
		geo.tiles <- tile.100[[1]]
		geo.meta <- lapply(tile.100[[2]], function(x) x[])

		# Check number of records for each tiles
		gbif.tiles <-
		vapply(
			geo.tiles,
			function(x) {
				gt.out <- rgbif::occ_count(
					taxonKey = sp.key,
					hasCoordinate = !no_xy,
					hasGeospatialIssue = FALSE,
					geometry = x
				)
				return(gt.out)
			},
			FUN.VALUE = numeric(1)
		)

		# We make sure that each tile will be fragmented enough
		# (i.e., < 10'000 species records each)
		keep.tiles <-
		lapply(seq_along(gbif.tiles), function(x)
		{
			# Which tile
			tile.n <- gbif.tiles[x]

			# Return NULL or POLYGON if 0 records or > 10'000 records
			if (tile.n == 0) {return(NULL)}
			if (tile.n < 10000) {return(geo.tiles[x])}

			# Keep fragmenting if > 10'000 records
			new.geo <- list(geo.meta[[x]])
			pol.store <- list()
			while (length(new.geo) != 0)
			{
				m.process <-
				lapply(seq_along(new.geo), function(y){

          Sys.sleep(3)

					# Convert to extent and create 5 new micro tiles
					new.ext <- terra::ext(new.geo[[y]])
					micro.100 <- make_tiles(new.ext, ntiles = 5, sext = TRUE)

					# If micro.100 is NULL return NULL (in case of an incorrect extent)
					if (is.null(micro.100)) {
						return(list(NULL,NULL))
					}

					# Continue
					micro.tiles <- micro.100[[1]]
					micro.meta <- lapply(micro.100[[2]], function(z) z[])

					# And count records again
					gbif.micro <-
						vapply(
							micro.tiles,
							function(z) {
								gt.out <- try(
									rgbif::occ_count(
										taxonKey = sp.key,
										hasCoordinate = !no_xy,
										hasGeospatialIssue = FALSE,
										geometry = z
									),
								silent = TRUE
							)
						return(gt.out)
						},
						FUN.VALUE = numeric(1)
					)

					# Transform in 0 when the polygon is too small for occ_search
					gbif.micro[
						vapply(
							gbif.micro,
							function(z) grepl("Error", z),
							FUN.VALUE = logical(1)
						)
					] <- 0
					gbif.micro <- as.numeric(gbif.micro)

					# Store if tiles with < 10'000 observations found
					goody <- micro.tiles[gbif.micro < 10000 & gbif.micro != 0]
					bady <- micro.meta[!gbif.micro < 10000 & gbif.micro != 0]

          # Return
					return(list(goody,bady))
				})

				# Restructure
				all.good <- lapply(m.process,function(y) y[[1]])
				all.good[vapply(all.good, is.null, logical(1))] <- NULL
				all.bad <- lapply(m.process, function(y) y[[2]])
				all.good[vapply(all.bad, is.null, logical(1))] <- NULL

				# Store if no blocks with > 10'000 observations found
				pol.store <- c(pol.store, unlist(all.good, recursive = FALSE))
				new.geo <- unlist(all.bad,recursive = FALSE)

				# Finally remove too small windows (for "invisible" duplicated records)
				inf.id <- vapply(
					new.geo,
					function(y) (y[2] - y[1]) < 1e-7,
					FUN.VALUE = logical(1)
				)
				new.geo <- new.geo[!inf.id]
			}
			return(unlist(pol.store, recursive = FALSE))
		})

		# Final list of POLYGONS with n observations < 10'000 observations
		geo.ref <- unlist(keep.tiles)
	}
	

	######################################################
	#################### API Search ######################
	######################################################


	# Information messages
	if (occ_samp != 10000) {
		cat("...GBIF records of",sp_name,": download of sample starting...","\n")
	} else {
		cat("...GBIF records of",sp_name,": download starting...","\n")
	}

	# Run the gbif search with the acceptedName per chosen tiles
	batch.search <-
	lapply(seq_along(geo.ref), function(x)
	{
    Sys.sleep(3)

		## Try the download first: may be request overload problems
		go.tile <- geo.ref[x]
		gbif.search <- if (should_use_occ_download) {
			
			## Try to use parameter creds if provided, otherwise use env variables
			user <- if (is.null(occ_download_user)) {
				Sys.getenv("GBIF_USER")
			} else {
				occ_download_user
			}
			pwd <- if (is.null(occ_download_pwd)) {
				Sys.getenv("GBIF_PWD")
			} else {
				occ_download_pwd
			}
			email <- if (is.null(occ_download_email)) {
				Sys.getenv("GBIF_EMAIL")
			} else {
				occ_download_email
			}

			if (any(c(user, pwd, email) %in% "")) {
				stop(paste(
					"GBIF credentials missing...Please provide the parameters",
					"'occ_download_user','occ_download_pwd' and 'occ_download_email'",
					"or set as environment variables GBIF_USER, GBIF_PWD",
					"and GBIF_EMAIL..."
				))
			}

			req_id <- rgbif::occ_download(
				rgbif::pred("taxonKey", sp.key),
				rgbif::pred("hasCoordinate", !no_xy),
				rgbif::pred("hasGeospatialIssue", FALSE),
				rgbif::pred_within(go.tile),
				format = "SIMPLE_CSV",
				curlopts = list(http_version=2),
				user = user,
				pwd = pwd,
				email = email
			)
			
			cat(">>> Download request ID:", req_id, "\n")
			
			rgbif::occ_download_wait(req_id,
				status_ping = 5,
				curlopts = list(http_version=2),
				quiet = FALSE)
			download <- rgbif::occ_download_get(req_id)
			try(rgbif::occ_download_import(download), silent = FALSE)
		
		} else {
    	cat("\r", "----------------- #", 
    		x, " (", round(x * 100/length(geo.ref), 2),
    		"%...)\033[K", sep="")
			try(
				rgbif::occ_search(taxonKey = sp.key,
								limit = occ_samp,
								hasCoordinate = !no_xy,
								hasGeospatialIssue = FALSE,
								geometry = go.tile), silent = TRUE)
		}

		# If problems, just rerun with while with n attempts, otherwise return NULL
		if (!should_use_occ_download && class(gbif.search) %in% "try-error") {
			print(gbif.search[1])
			warning(
				"\n",
				"GBIF query overload or rgbif package error [taxonKey=",
				sp.key,"]...","\n",sep=""
			)

			# While
			if (class(gbif.search) %in% "try-error") {
				j <- 0
				while (class(gbif.search) %in% "try-error" & j < ntries)
				{
					cat("\n","Attempt", j + 1, "...", "\n")
					j <- j + 1
					gbif.search <- try(
						rgbif::occ_search(taxonKey = sp.key,
										limit = occ_samp,
										hasCoordinate = !no_xy,
										hasGeospatialIssue = FALSE,
										geometry = go.tile), silent = TRUE)
					  Sys.sleep(3)
				}

				if (class(gbif.search) %in% "try-error") {
					if (error_skip){
						cat("\n","Attempts to download failed...Returning no results")
						return(e.output)

					} else {
						stop(
							"\n","ERROR (not skipped) for [taxonKey=",sp.key,"]...",
							"\n",sep=""
						)
					}
				}
			}
		}

		if (should_use_occ_download) gbif.search$data <- gbif.search

		# If no results
		if (is.null(gbif.search$data)){
			return(e.output)
		
		} else {
			# Convert to a data.frame is needed
			if (class(gbif.search$data)[1] != "data.frame") {
				gbif.search <- as.data.frame(gbif.search$data)
			}

			# Reordering data.frame and correcting if missing columns
			gbif.reorder <- gbif.search[, order(names(gbif.search))]
			missing.col <- gbif.info[!(gbif.info %in% names(gbif.reorder))]
			gbif.reorder[, missing.col] <- NA
			gbif.correct <- gbif.reorder[, order(names(gbif.reorder))]

			# Remove columns we don't want
			return(gbif.correct[, gbif.info])
		}
	})

	# Combine all results in one data.frame
	gbif.compile <- do.call("rbind",batch.search)

	# Keep specific fields (just for safety)
	gbif.compile <- gbif.compile[, gbif.info]


	#############################################################
	#################### Records filtering ######################
	#############################################################


	#### 5) Select records according to year range
	cat("---> Time period selection...","\n")

	id.year <- gbif.correct$year >= min(time_period) &
								gbif.correct$year <= max(time_period)
	id.year[is.na(id.year)] <- TRUE
	gbif.correct <- gbif.correct[id.year, ]

		# Removal summary
	if (any(names(table(id.year)) %in% FALSE)){
		removed <- table(id.year)[1]
		cat("Records removed:",removed,"\n")

	} else {
		cat("Records removed:",0,"\n")
	}

	#### 6) Remove records with identical xy
	if (nrow(gbif.correct) > 0){
		if (!identic_xy){
		
			cat("---> Removal of identical xy records...","\n")

			id.diff <- c(!(abs(gbif.correct[, "decimalLatitude"]) ==
				abs(gbif.correct[, "decimalLongitude"])))
			gbif.correct <- gbif.correct[id.diff, ]

			# Removal summary
			if (any(names(table(id.diff)) %in% FALSE)){
				removed <- table(id.diff)[1]
				cat("Records removed:",removed,"\n")

			} else {
				cat("Records removed:",0,"\n")
			}
		}
	}

	#### 7) Remove wrongly lon/lat converted (only if > 50 records in one dataset)
	gbif.correct <- as.data.frame(gbif.correct)
	gbif.nrow <- nrow(gbif.correct)

	if (nrow(gbif.correct) > 0){
		if (!wConverted_xy){

			cat("---> Removal of wrong lon/lat converted records...","\n")

			# Keep dataset with NAs
			gbif.na <- gbif.correct[is.na(gbif.correct$datasetKey),]

			# Summary of datasets & cd_ddmm parameter choice
			gbif.datasets <- names(table(gbif.correct$datasetKey))
			if (nrow(gbif.correct) < 10000) {
				mat.size <- 100
			} else if (nrow(gbif.correct) < 1e6) {
				mat.size <- 1000
			} else {
				mat.size <- 10000
			}

			# Apply correction if the dataset > 50 records
			gbif.ddmm <-
			lapply(seq_along(gbif.datasets), function(x)
			{
				gbif.dataset <- gbif.correct[gbif.correct$datasetKey %in% gbif.datasets[x],]
				if (nrow(gbif.dataset) > 50){
					diff.chosen <-
						diff.records[[2]][which(nrow(gbif.dataset) < diff.records[[1]])[1]]
					gbif.dataset <- suppressWarnings(
						CoordinateCleaner::cd_ddmm(
							x = gbif.dataset,
							lon = "decimalLongitude",
							lat = "decimalLatitude",
							ds = "datasetKey",
							mat_size = mat.size,
							diff = diff.chosen,
							min_span = xy.span
						)
					)
				}
				return(gbif.dataset)
			})

			gbif.Wna <- do.call("rbind", gbif.ddmm)
			gbif.correct <- rbind(gbif.na, gbif.Wna)

			# Removal summary
			cat("Records removed:",gbif.nrow - nrow(gbif.correct),"\n")
		}
	}

	#### 8) Remove raster centroids (only if > 100 records in one dataset)
	gbif.nrow <- nrow(gbif.correct)

	if (nrow(gbif.correct) > 0){
		if (!centroids){

			cat("---> Removal of raster centroids...","\n")

			# Keep dataset with NAs
			gbif.na <- gbif.correct[is.na(gbif.correct$datasetKey), ]
			
			# Summary of datasets & cd_ddmm parameter choice
			gbif.datasets <- names(table(gbif.correct$datasetKey))

			# Apply correction if the dataset > 100 records
			gbif.round <-
			lapply(seq_along(gbif.datasets), function(x){

				gbif.dataset <- gbif.correct[gbif.correct$datasetKey %in% gbif.datasets[x],]
				if (nrow(gbif.dataset) > 100){
					gbif.temp <- suppressWarnings(
						try(CoordinateCleaner::cd_round(
							x = gbif.dataset,
							lon = "decimalLongitude",
							lat = "decimalLatitude",
							ds = "datasetKey",
							graphs = FALSE,
							...
							),
						silent = TRUE)
					)
					
					if (class(gbif.temp) %in% "try-error") {
						return(gbif.dataset)

					} else {
						gbif.dataset <- gbif.temp
					}
				}
				return(gbif.dataset)
			})

			gbif.Wna <- do.call("rbind", gbif.round)
			gbif.correct <- rbind(gbif.na, gbif.Wna)

			# Removal summary
			cat("Records removed:",gbif.nrow - nrow(gbif.correct),"\n")
		}
	}
	if (nrow(gbif.correct) == 0) {
		cat("No records left after selection...","\n")
		return(e.output)

	} else {
		proper.output <- getGBIF(cbind(input_search = sp_name,gbif.correct))
		return(proper.output)
	}
}
