### =========================================================================
### get_gbif
### =========================================================================
#' Download and Filter GBIF Occurrences for Spatial Analyses
#'
#' Retrieve GBIF occurrence records for a focal taxon, harmonize the supplied
#' name against the GBIF backbone taxonomy, optionally split large downloads
#' into smaller tiles, and apply a series of post-download filters designed for
#' spatial analyses and range mapping.
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
#' @param has_xy Logical. If \code{TRUE} (default), keep only records with
#' coordinates. If \code{FALSE}, keep only records without coordinates. If
#' \code{NULL}, keep all records.
#' @param spatial_issue Logical. If \code{FALSE} (default), keep only records
#' without geospatial issues. If \code{TRUE}, keep only records with
#' geospatial issues. If \code{NULL}, keep all records.
#' @param grain Numeric study grain in kilometers. Default is \code{100}. The
#' value is used to filter records by \code{coordinateUncertaintyInMeters} and
#' by the number of reported coordinate decimals.
#' @param duplicates Logical. Should duplicate records be kept? Default is
#' \code{FALSE}.
#' @param absences Logical. Should absence records be kept? Default is
#' \code{FALSE}.
#' @param basis Character vector giving the accepted bases of record. The
#' default keeps commonly used occurrence-oriented record types and excludes
#' specimen- and unknown-based records.
#' @param establishment Character vector giving accepted
#' \code{degreeOfEstablishment} values. The default keeps native and broadly
#' established records, together with records lacking that field.
#' @param add_infos Optional character vector of additional GBIF occurrence
#' fields to append to the default output.
#' @param time_period Numeric vector of length two giving the year range to
#' retain. Default is \code{c(1000, 3000)}. Records with missing years are
#' kept.
#' @param identic_xy Logical. Should records with identical longitude-latitude
#' pairs be kept? Default is \code{FALSE}.
#' @param wConverted_xy Logical. Should records that appear to be incorrectly
#' converted from degree-minute notation be kept? Default is \code{TRUE}. If
#' \code{FALSE}, an approximate version of
#' \code{CoordinateCleaner::cd_ddmm()} is applied.
#' @param centroids Logical. Should records located on raster centroids be
#' kept? Default is \code{FALSE}. If \code{FALSE},
#' \code{CoordinateCleaner::cd_round()} is used.
#' @param ntries Numeric number of download attempts before giving up after
#' GBIF or \code{rgbif} errors. Default is \code{10}.
#' @param error_skip Logical. If \code{TRUE}, return an empty result when all
#' download attempts fail.
#' @param occ_samp Numeric maximum number of GBIF occurrences sampled per
#' geographic tile. Default is \code{10000}.
#' @param should_use_occ_download Logical. If \code{TRUE}, use
#' \code{rgbif::occ_download()} instead of \code{rgbif::occ_search()}.
#' This requires GBIF credentials. Default is \code{FALSE}.
#' @param occ_download_user Character GBIF username. Required if
#' \code{should_use_occ_download = TRUE}.
#' @param occ_download_pwd Character GBIF password. Required if
#' \code{should_use_occ_download = TRUE}.
#' @param occ_download_email Character GBIF email address. Required if
#' \code{should_use_occ_download = TRUE}.
#' @param verbose Logical. Should progress messages be printed? Default is
#' \code{TRUE}.
#' @param ... Additional arguments passed to
#' \code{CoordinateCleaner::cd_round()}.
#' @details The function follows the same taxonomic matching logic used by the
#' GBIF website and retrieves records linked to the accepted name and its
#' synonyms. Internally, the input name is first resolved against the GBIF
#' backbone taxonomy. If the matched name is a synonym, the download is carried
#' out with the corresponding accepted GBIF taxon key rather than the synonym
#' key itself.
#'
#' This means that \code{get_gbif()} harmonizes the taxon concept used for the
#' query, but it does not rewrite every returned record to a single final name.
#' The output keeps the original GBIF taxonomic fields for each record,
#' including \code{scientificName}, while also retaining the harmonized fields
#' \code{acceptedScientificName}, \code{acceptedTaxonKey}, and
#' \code{taxonomicStatus}.
#' 
#' If the requested extent contains many records, the query is fragmented into
#' a set of tiles so that individual API calls remain manageable. By default,
#' tiles are refined until they contain fewer than about 10,000 records, which
#' is faster and more robust than relying on a single large request.
#'
#' Post-download filtering can remove records outside the study extent,
#' duplicated coordinates, absences, unwanted bases of record, records with low
#' spatial precision, suspicious coordinate conversions, and raster-centroid
#' records. The \code{grain} argument controls both coordinate-uncertainty
#' filtering and the minimum number of coordinate decimals that must be
#' reported.
#'
#' Decimal filtering follows these thresholds:\cr
#' - if 110 km > \code{grain} \eqn{\ge}{>=} 11 km, coordinates with at least
#' 1 decimal are kept\cr
#' - if 11 km > \code{grain} \eqn{\ge}{>=} 1100 m, coordinates with at least
#' 2 decimals are kept\cr
#' - if 1100 m > \code{grain} \eqn{\ge}{>=} 110 m, coordinates with at least
#' 3 decimals are kept\cr
#' - if 110 m > \code{grain} \eqn{\ge}{>=} 11 m, coordinates with at least
#' 4 decimals are kept\cr
#' - if 11 m > \code{grain} \eqn{\ge}{>=} 1.1 m, coordinates with at least
#' 5 decimals are kept
#' @return A \code{getGBIF} object (a data.frame subclass) containing the
#' requested occurrence data. The object may also carry the attributes
#' \code{filter_log}, which records how many rows were removed at each
#' filtering step, and \code{no_xy}, which stores records without coordinates
#' when they are requested.
#'
#' Taxonomic harmonization is reflected in the output columns rather than by
#' replacing all names with a single accepted label. In particular,
#' \code{scientificName} stores the record-level GBIF name, whereas
#' \code{acceptedScientificName}, \code{acceptedTaxonKey}, and
#' \code{taxonomicStatus} expose the accepted GBIF interpretation used by the
#' query.
#' @references
#' Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger,
#' D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and
#' land cover on plant species distribution in the European Alps. Ecological
#' monographs, 91(2), e01433. 10.1002/ecm.1433
#'
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the
#' global biodiversity information facility API. 10.5281/zenodo.6023735
#'
#' Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C.,
#' Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized
#' cleaning of occurrence records from biological collection databases.
#' Methods in Ecology and Evolution, 10(5), 744-751. 10.1111/2041-210X.13152
#'
#' Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7."
#' (2022). Terra - CRAN
#' @seealso \code{get_status()} to inspect the accepted name and synonym mapping
#' returned by GBIF; \code{rgbif} for more general GBIF retrieval workflows;
#' and \code{CoordinateCleaner} for more extensive occurrence cleaning.
#' @example inst/examples/get_gbif_help.R
#' @importFrom terra ext vect
#' @importFrom rgbif name_backbone occ_search
#' @importFrom CoordinateCleaner cd_ddmm cd_round
#' @importFrom stats na.omit
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
					has_xy = TRUE,
					spatial_issue = FALSE,
					grain = 100,
					duplicates = FALSE,
					absences = FALSE,
					basis =  c('OBSERVATION', 'HUMAN_OBSERVATION', 'MACHINE_OBSERVATION',
						'OCCURRENCE', 'MATERIAL_CITATION', 'MATERIAL_SAMPLE','LITERATURE'),
					establishment = c('native','casual','released','reproducing',
						'established','colonising','invasive','widespreadInvasive'),
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
					verbose = TRUE,
					...) {
	
	######################################################
	### Stop messages
	######################################################


	# General
	check_character_vector(sp_name, "sp_name")
	check_logical(search, "search")
	check_numeric(conf_match, "conf_match")
	check_numeric(grain, "grain")
	check_logical(duplicates, "duplicates")
	check_logical(absences, "absences")
	check_character_vector(establishment, "establishment")
	check_numeric_range(time_period, "time_period", 2)
	check_logical(identic_xy, "identic_xy")
	check_logical(wConverted_xy, "wConverted_xy")
	check_logical(centroids, "centroids")
	check_numeric(ntries, "ntries")
	check_logical(error_skip, "error_skip")
	check_numeric(occ_samp, "occ_samp")
	check_logical(should_use_occ_download, "should_use_occ_download")

	# Basis of records
	allb <- c('OBSERVATION', 'HUMAN_OBSERVATION', 'MACHINE_OBSERVATION',
		'MATERIAL_CITATION', 'MATERIAL_SAMPLE', 'PRESERVED_SPECIMEN',
		'FOSSIL_SPECIMEN', 'LIVING_SPECIMEN', 'LITERATURE', 'UNKNOWN', 'OCCURRENCE')
	if (!all(basis %in% allb) || is.null(basis)) {
		stop("Wrong given 'basis', see function help --> ?get_gbif()")
	}

	# Establishment
	alle <- c('native', 'casual', 'released', 'reproducing', 'established',
		'colonising', 'invasive', 'widespreadInvasive', 'managed', 'captive',
		'cultivated', 'released', 'unestablished', 'failing')
	if (!all(establishment %in% alle) || is.null(establishment)) {
		stop("Wrong given 'establishment', see function help --> ?get_gbif()")
	}

	# Add_infos
	if (all(!is.null(add_infos))) {
		if (all(is.na(add_infos))) {
			stop("NA values not allowed for given 'add_infos'")
		}
	}

	# Geo class
	spatial.class <- c("Extent", "SpatExtent", "SpatialPolygon",
		"SpatialPolygonDataframe", "SpatVector", "sf")
	if (!is.null(geo)) {
		if (!class(geo)[1] %in% spatial.class) {
			stop("Wrong 'geo' class (not a spatial object)...")
		}
	}


	######################################################
	### Parameter config
	######################################################


	# For precision and 'cd_ddmm'
		# Convert grain to needed units
	grain_km  <- grain
	grain_m   <- grain_km * 1000

	# For fields
	gbif.info <- c('taxonKey', 'scientificName', 'acceptedTaxonKey',
		'acceptedScientificName', 'individualCount', 'occurrenceStatus',
		'establishmentMeans', 'degreeOfEstablishment', 'decimalLatitude',
		'decimalLongitude', 'basisOfRecord', 'coordinateUncertaintyInMeters',
		'countryCode', 'country', 'year', 'datasetKey', 'institutionCode',
		'publishingOrgKey', 'taxonomicStatus', 'taxonRank', add_infos)
	gbif.info <- gbif.info[order(gbif.info)]

	# Create an empty ouptut
	e.output <- getGBIF(
		data.frame(matrix(
			ncol = length(gbif.info),
			nrow = 0,
			dimnames = list(NULL, gbif.info))
		)
	)

	# Summary helper
	summary_log <- data.frame(
		step = character(),
		removed = integer(),
		remaining = integer(),
		stringsAsFactors = FALSE
	)

	# Print helper
	vcat <- function(..., sep = "") {
		if (isTRUE(verbose)) cat(..., sep = sep)
	}


	######################################################
	### Backbone harmonization and N records
	######################################################


	if (!search) {
	# Search input name via fuzzy match and direct search
	bsearch <- rgbif::name_backbone(
		sp_name,
		rank = rank,
		phylum = phylum,
		class = class,
		order = order,
		family = family,
		verbose = FALSE,
		strict = FALSE
	)

	} else {
		# Search input name via strict match and refined search
		bsearch <- rgbif::name_backbone(
			sp_name,
			verbose = TRUE,
			strict = TRUE
		)

		q.crit <- !vapply(
			list(rank, phylum, class, order, family),
			is.null,
			logical(1)
		)

		# Filter by the supplied criteria if results are available
		if (!bsearch$matchType[1] %in% "NONE") { 
			if (any(q.crit)) {
				id.crit <- c("rank", "phylum", "class", "order", "family")[q.crit]
				p.crit <- unlist(list(rank, phylum, class, order, family)[q.crit])
				n.test <- id.crit %in% names(bsearch)

				if (any(n.test)) {
					# Selecting which
					id.crit2 <- id.crit[n.test]
					p.crit2 <- p.crit[n.test]

					# Apply the requested criteria
					for (i in seq_along(id.crit2)) {
						bsearch <- bsearch[c(bsearch[, id.crit2[i]])[[1]] %in% p.crit2[i],]

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
		if (nrow(bsearch) > 1) {
			if (all(!bsearch$rank %in% c("SPECIES", "SUBSPECIES", "VARIETY"))) {
				cat("Not match found...", "\n")
				return(e.output)

			} else {
				s.keep <- bsearch[bsearch$rank %in%
					c("SPECIES" ,"SUBSPECIES" ,"VARIETY"),]
				s.keep <- s.keep[s.keep$status %in% c("ACCEPTED", "SYNONYM"),]
				if (nrow(s.keep) == 0) {
					cat("Not match found...", "\n")
					return(e.output)

				} else if (nrow(s.keep) > 1) {
					# If only subspecies and variety are found, prioritize one
					if (all(s.keep$rank %in% c("VARIETY", "SUBSPECIES"))) {
						if ("var." %in% strsplit(sp_name," ")[[1]]) {
							bsearch <- s.keep[s.keep$rank %in% "VARIETY", ]
							if (nrow(bsearch) == 0) {
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
						length(unique(key.test[, 1])) == 1) {
						bsearch <- bsearch[bsearch$status %in% "ACCEPTED", ]
					}

				} else {
					bsearch <- s.keep
				}

				# If not the same species overall return empty
				s.usp <- length(unique(bsearch$speciesKey)) == 1
				if (!s.usp) {
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
	if (bsearch$status %in% "SYNONYM") {
		sp.key <- bsearch$acceptedUsageKey

	} else {
		sp.key <- bsearch$usageKey
	}


	#############################################################
	### Convert to POLYGON
	#############################################################


	if (!is.null(geo)) {
		if (inherits(geo, "sf")) geo <- terra::vect(geo)
		if (!inherits(geo, "SpatExtent")) geo <- terra::ext(geo)

		geo.ref <- list(
			paste0(
				"POLYGON((", geo$xmin, " ", geo$ymin, ", ",
					geo$xmax, " ", geo$ymin, ", ",
					geo$xmax, " ", geo$ymax, ", ",
					geo$xmin, " ", geo$ymax, ", ",
					geo$xmin, " ", geo$ymin, "))"
			)
		)
	
	} else {
		geo <- terra::ext()
		geo.ref <- list(NULL)
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
		geometry = geo.ref[[1]]
	)

	# Print summary
	l <- sprintf(
		"%-29s : %10d",
		c("Total number (all records)", "Kept records"),
		c(gbif.total, gbif.records)
	)
	w <- max(nchar(l)) + 4

	if (verbose) {
		cat(
			"|",strrep("-",w-2),"|\n| ",
			paste(l,collapse=" |\n| "),
			" |\n|",strrep("-",w-2),
			"|\n",sep=""
		)

		fmt <- function(x) if (is.null(x)) "NULL" else as.character(x)
		cat("| Kept records according to parameters:\n")

		# Print additional information depending if global or regional
		if (is.null(geo.ref[[1]])) {
			cat(sprintf(
					"| spatial_issue = %s, has_xy = %s",
					fmt(spatial_issue),
					fmt(has_xy)
				)
			)
		} else {
			cat(sprintf(
					paste0(
						"| spatial_issue = %s, has_xy = TRUE by default ",
						"('geo' was set)"
					),
					fmt(spatial_issue)
				)
			)
		}
	}

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

		vcat("\n...(> 10'000 records) retrieving tiles...\n")

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
					hasCoordinate = has_xy,
					hasGeospatialIssue = spatial_issue,
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
										hasCoordinate = has_xy,
										hasGeospatialIssue = spatial_issue,
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
		vcat("...GBIF records of ", sp_name, ": download of sample starting...\n")
	} else {
		vcat("...GBIF records of ", sp_name, ": download starting...\n")
	}

	# Run the gbif search with the acceptedName per chosen tiles
	batch.search <-
	lapply(seq_along(geo.ref), function(x)
	{
		Sys.sleep(3)

		## Try the download first: may be request overload problems
		go.tile <- geo.ref[[x]]
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
				rgbif::pred("hasCoordinate", has_xy),
				rgbif::pred("hasGeospatialIssue", spatial_issue),
				rgbif::pred_within(go.tile),
				format = "SIMPLE_CSV",
				curlopts = list(http_version=2),
				user = user,
				pwd = pwd,
				email = email
			)
			
			vcat(">>> Download request ID:", req_id, "\n")
			
			rgbif::occ_download_wait(
				req_id,
				status_ping = 5,
				curlopts = list(http_version=2),
				quiet = FALSE
			)
			download <- rgbif::occ_download_get(req_id)
			try(rgbif::occ_download_import(download), silent = FALSE)
		
		} else {
			vcat(
				"\r", "------------- #", 
				x, " (", round(x * 100/length(geo.ref), 2),
				"%...)\033[K",
				sep=""
			)
			try(
				rgbif::occ_search(
					taxonKey = sp.key,
					limit = occ_samp,
					hasCoordinate = has_xy,
					hasGeospatialIssue = spatial_issue,
					geometry = go.tile
				),
				silent = TRUE
			)
		}

		# If problems, just rerun with while with n attempts, otherwise return NULL
		if (!should_use_occ_download && class(gbif.search) %in% "try-error") {
			print(gbif.search[1])
			warning(
				"\n","GBIF query overload or rgbif package error [taxonKey=",
				sp.key,"]...","\n",sep=""
			)

			# While
			if (class(gbif.search) %in% "try-error") {
				j <- 0
				while (class(gbif.search) %in% "try-error" & j < ntries)
				{
					vcat("Attempt", j + 1, "...", "\n")
					j <- j + 1
					gbif.search <- try(
						rgbif::occ_search(
							taxonKey = sp.key,
							limit = occ_samp,
							hasCoordinate = has_xy,
							hasGeospatialIssue = spatial_issue,
							geometry = go.tile
						),
						silent = TRUE
					)
					Sys.sleep(3)
				}

				if (class(gbif.search) %in% "try-error") {
					if (error_skip){
						cat("Attempts to download failed...Returning no results","\n")
						return(e.output)

					} else {
						stop(
							"ERROR (not skipped) for [taxonKey=",sp.key,"]...",
							"\n",sep=""
						)
					}
				}
			}
		}

		if (should_use_occ_download) gbif.search$data <- gbif.search

		# If no results
		if (is.null(gbif.search$data)) {
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

	# Initial number
	n_initial <- nrow(gbif.compile)

	# Has xy / not
	gbif.no_xy <- gbif.compile[
		is.na(gbif.compile$decimalLatitude) |
		is.na(gbif.compile$decimalLongitude),
	]
	gbif.xy <- gbif.compile[
		!is.na(gbif.compile$decimalLatitude) &
		!is.na(gbif.compile$decimalLongitude),
	]

	#### 1) Grain filtering
		
	# GBIF uncertainty
	id.certain <- gbif.xy$coordinateUncertaintyInMeters <= grain_m / 2
	id.certain[is.na(id.certain)] <- TRUE
	gbif.correct <- gbif.xy[id.certain, ]

	# GBIF lon/lat decimals
	if (grain_km < 110) {

		# Remove latitude or/and longitude with no decimals
		lonlat.format <- data.frame(
			lat = as.character(gbif.correct$decimalLatitude),
			lon = as.character(gbif.correct$decimalLongitude)
		)

		# Require decimals on BOTH lat and lon
		id.deci <- grepl("\\.", lonlat.format$lat) & grepl("\\.", lonlat.format$lon)
		gbif.correct <- gbif.correct[id.deci, ]

		grain_deg <- grain_m / 111320 # 111,320 meters ≈ 1 degree
		deci.chosen <- max(0, ceiling(-log10(grain_deg)))

		# Keep coordinates compatible with the input 'grain'
		declat <- gsub(".*\\.", "", lonlat.format[id.deci,1])
		declon <- gsub(".*\\.", "", lonlat.format[id.deci,2])
		id.grain <- nchar(declon) >= deci.chosen & nchar(declat) >= deci.chosen
		gbif.correct <- gbif.correct[id.grain, ]
	}

	# Removal summary
	n_after <- nrow(gbif.correct)
	summary_log <- log_step(
		log = summary_log, 
		step_name = "Grain filtering", 
		before = n_initial,
		after = n_after
	)


	#### 2) Removing xy duplicates
	if (!duplicates) {

		n_before <- nrow(gbif.correct)
		
		id.dup <- !duplicated(gbif.correct[, c("decimalLongitude","decimalLatitude")])
		gbif.correct <- gbif.correct[id.dup, ]

		# Removal summary
		n_after <- nrow(gbif.correct)
		summary_log <- log_step(
			log = summary_log, 
			step_name = "Duplicated records", 
			before = n_before,
			after = n_after
		)
	}


	#### 3) Removing absences
	if (!absences) {

		n_before <- nrow(gbif.correct)
		
		id.abs <- !(gbif.correct$individualCount %in% 0 |
					gbif.correct$occurrenceStatus %in% "ABSENT")
		gbif.correct <- gbif.correct[id.abs, ]

		# Removal summary
		n_after <- nrow(gbif.correct)
		summary_log <- log_step(
			log = summary_log, 
			step_name = "Absence records", 
			before = n_before,
			after = n_after
		)
	}


	#### 4) Select basis of records
	n_before <- nrow(gbif.correct)
	id.basis <- gbif.correct$basisOfRecord %in% basis
	gbif.correct <- gbif.correct[id.basis, ]
	
	# Removal summary
	n_after <- nrow(gbif.correct)
	summary_log <- log_step(
		log = summary_log, 
		step_name = "Basis selection", 
		before = n_before,
		after = n_after
	)


	#### 5) Select establishment of records
	n_before <- nrow(gbif.correct)
	id.esta <- gbif.correct$degreeOfEstablishment %in% establishment |
					is.na(gbif.correct$degreeOfEstablishment)
	gbif.correct <- gbif.correct[id.esta, ]
	
	# Removal summary
	n_after <- nrow(gbif.correct)
	summary_log <- log_step(
		log = summary_log, 
		step_name = "Establishment selection", 
		before = n_before,
		after = n_after
	)


	#### 6) Select records according to year range
	n_before <- nrow(gbif.correct)
	id.year <- gbif.correct$year >= min(time_period) &
					gbif.correct$year <= max(time_period)
	id.year[is.na(id.year)] <- TRUE
	gbif.correct <- gbif.correct[id.year, ]
	
	# Removal summary
	n_after <- nrow(gbif.correct)
	summary_log <- log_step(
		log = summary_log, 
		step_name = "Time frame", 
		before = n_before,
		after = n_after
	)

	
	#### 7) Remove records with identical xy
	if (nrow(gbif.correct) > 0) {
		if (!identic_xy){
			n_before <- nrow(gbif.correct)

			id.diff <- !(abs(gbif.correct[, "decimalLatitude"]) ==
						 abs(gbif.correct[, "decimalLongitude"]))
			gbif.correct <- gbif.correct[id.diff, ]

			# Removal summary
			n_after <- nrow(gbif.correct)
			summary_log <- log_step(
				log = summary_log, 
				step_name = "Identical records", 
				before = n_before,
				after = n_after
			)
		}
	}


	#### 8) Remove lon/lat missconverted
	if (!wConverted_xy && nrow(gbif.correct) > 0) {
		n_before <- nrow(gbif.correct)
		gbif.correct <- as.data.frame(gbif.correct)

		# Keep records with NA datasetKey
		gbif.na <- gbif.correct[is.na(gbif.correct$datasetKey), ]

		# Unique dataset keys (excluding NA)
		gbif.datasets <- unique(na.omit(gbif.correct$datasetKey))

		# Apply correction only if dataset > 50 records
		gbif.ddmm <-
		lapply(gbif.datasets, function(ds) {
			gbif.dataset <- gbif.correct[gbif.correct$datasetKey == ds, ]
			if (nrow(gbif.dataset) > 50) {
				# Set
				n <- nrow(gbif.dataset)
				# Matrix size
				mat.size <- if (n < 10000) {
					100
				} else if (n < 1e6) {
					1000
				} else {
					10000
				}
				# Dataset spatial extent
				lon_span <- diff(range(gbif.dataset$decimalLongitude, na.rm = TRUE))
				lat_span <- diff(range(gbif.dataset$decimalLatitude, na.rm = TRUE))
				dataset_span <- max(lon_span, lat_span)
				# Minimum span
				min.span <- max(1, min(5, dataset_span * 0.1))
				# Decimal imbalance threshold
				diff.param <- 1
				# Function
				gbif.dataset <- suppressWarnings(
					CoordinateCleaner::cd_ddmm(
						x = gbif.dataset,
						lon = "decimalLongitude",
						lat = "decimalLatitude",
						ds = "datasetKey",
						mat_size = mat.size,
						diff = diff.param,
						min_span = min.span
					)
				)
			}
			return(gbif.dataset)
		})
		gbif.Wna <- if (length(gbif.ddmm) > 0) {
			do.call("rbind", gbif.ddmm)
		} else {
			NULL
		}
		gbif.correct <- if (!is.null(gbif.Wna)) {
			rbind(gbif.na, gbif.Wna)
		} else {
			gbif.na
		}
		
		# Removal summary
		n_after <- nrow(gbif.correct)
		summary_log <- log_step(
			log = summary_log, 
			step_name = "Missconverted records", 
			before = n_before,
			after = n_after
		)
	}


	#### 9) Remove raster centroids
	if (!centroids && nrow(gbif.correct) > 0) {
		n_before <- nrow(gbif.correct)
		gbif.correct <- as.data.frame(gbif.correct)

		# Keep records with NA datasetKey
		gbif.na <- gbif.correct[is.na(gbif.correct$datasetKey), ]

		# Unique non-NA datasets
		gbif.datasets <- unique(na.omit(gbif.correct$datasetKey))

		# Apply correction only if dataset > 100 records
		gbif.round <- lapply(gbif.datasets, function(ds) {
			gbif.dataset <- gbif.correct[gbif.correct$datasetKey == ds, ]

			if (nrow(gbif.dataset) > 100) {

				gbif.temp <- suppressWarnings(
					try(
						CoordinateCleaner::cd_round(
							x = gbif.dataset,
							lon = "decimalLongitude",
							lat = "decimalLatitude",
							ds = "datasetKey",
							graphs = FALSE,
							verbose = FALSE,
							...
						),
						silent = TRUE
					)
				)

				if (!inherits(gbif.temp, "try-error")) {
					gbif.dataset <- gbif.temp
				}
			}
			return(gbif.dataset)
		})	
		gbif.Wna <- if (length(gbif.round) > 0) {
			do.call("rbind", gbif.round)
		} else {
			NULL
		}
		gbif.correct <- if (!is.null(gbif.Wna)) {
			rbind(gbif.na, gbif.Wna)
		} else {
			gbif.na
		}

		# Removal summary
		n_after <- nrow(gbif.correct)
		summary_log <- log_step(
			log = summary_log, 
			step_name = "Raster centroids", 
			before = n_before,
			after = n_after
		)
	}

	if (nrow(gbif.correct) == 0) {
	    if (nrow(gbif.no_xy) > 0) {
	        if (verbose) {
	            cat("No spatial records left after filtering...\n")
	        }

	        proper.output <- getGBIF(
	            cbind(
	                input_search = sp_name,
	                gbif.no_xy
	            )
	        )

	        attr(proper.output, "filter_log") <- summary_log
	        attr(proper.output, "no_xy") <- gbif.no_xy
	        return(proper.output)

	    } else {
	        cat("No records left after filtering...\n")
	        return(e.output)
	    }
	}

	# ---- FINAL SUMMARY PRINT ----
	if (verbose) {
		if (nrow(summary_log) > 0) {
			max_step_width <- max(nchar(summary_log$step))
			max_num_width <- max(
				nchar(as.character(summary_log$removed)),
				nchar(as.character(summary_log$remaining))
			)
		} else {
			max_step_width <- 10
			max_num_width <- 10
		}
		
		box_width <- max_step_width + max_num_width + 20
		sep_dash  <- strrep("-", box_width)

		vcat("...Records (XY) filtering summary:\n")

		vcat(sep_dash, "\n")
		print(summary_log, row.names = FALSE)
		vcat("\n")
		vcat(sprintf("%-23s : %d\n", "Initial records", n_initial))
		vcat(sprintf("%-23s : %d\n", "Total removed", n_initial - nrow(gbif.correct)))
		vcat(sprintf("%-23s : %d\n", "Final records (XY)", nrow(gbif.correct)))
		vcat(sep_dash, "\n")
		vcat(sprintf("%-23s : %d\n", "Final records (no XY)", nrow(gbif.no_xy)))
	}

	# ---- RETURN AFTER PRINTING ----
	proper.output <- getGBIF(
		cbind(
			input_search = sp_name,
			gbif.correct
		)
	)

	attr(proper.output, "filter_log") <- summary_log
	attr(proper.output, "no_xy") <- gbif.no_xy
	return(proper.output)
}
