### =========================================================================
### get_gbif
### =========================================================================
#' Massively download and filter GBIF observations for sound spatial analyses
#'
#' Implement an user-friendly workflow to download and clean gbif taxa observations.
#' The function uses the rgbif R package but (1) implements the same search result 
#' found if www.gbif.org is employed i.e., based on the input taxa name, all species
#' records related to its accepted name and synonyms are extracted. The function
#' also (2) bypasses the rgbif hard limit on the number of records (100'000 max).
#' For this purpose, a dynamic moving window is created and used across the geographic
#' extent defined by the user. This window automatically fragments the specified
#' study area in succesive tiles of different sizes, until all tiles include < 100'000
#' observations. The function also (3) automatically applies a post-filtering of
#' observations based on the chosen resolution of the study/analysis and by partly
#' employing the CoordinateCleaner R package. Filtering options may be chosen and
#' involve several choices: study's extent, removal of duplicates, removal of absences,
#' basis of records selection, removal of invalid/uncertain xy coordinates (WGS84), time
#' period selection and removal of raster centroids. By default, the argument
#' hasGeospatialIssue in occ_data() (implemented rgbif function) is set to FALSE.
#' 
#' @param sp_name Character. Species name from which the user wants to retrieve all existing GBIF names
#' with associated taxonomy and IUCN status.
#' @param search Logical. If TRUE, the function will strictly look for the most relevant result, based
#' on a list of names given by rgbif (only species, subspecies and variety allowed here), and give an error
#' if name matching was impeded by synonym duplicates. If FALSE, the function will simply pick the first most
#' relevant name from the list (higher taxa level than species allowed here). Also, unlike search=TRUE,
#' fuzzy search (~approximative name match) is here allowed, and the 'rank', phylum', 'class', order'
#' and 'family' parameters are optionally used only if no convincing name match is found.
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
#' @param conf_match Numeric from 0 to 100. Determine the confidence threshold of match
#' of 'sp_name' with the GBIF backbone taxonomy. Default is 90.
#' @param geo Object of class Extent, SpatExtent, SpatialPolygon, SpatialPolygonDataframe,
#' or SpatVector (WGS84) to define the study's area extent. Default is NULL i.e. the whole globe.
#' @param grain Numeric. Specify in meters the study resolution. Used to
#' filter gbif records according to their (1) spatial uncertainties and (2) number of coordinate
#' decimals. Records with no information on coordinate uncertainties (column
#' 'coordinateUncertaintyInMeters') are be kept by default. See details.
#' @param duplicates Logical. Should duplicated records be kept?
#' @param absences Logical. Should absence records be kept?
#' @param no_xy Logical. Default is FALSE i.e., only records with coordinates are
#' downloaded. If TRUE, records with no coordinates are also downloaded.
#' @param basis Character. Which basis of records should be selected?
#' Available (old and new) are 'OBSERVATION', 'HUMAN_OBSERVATION', 'MACHINE_OBSERVATION',
#' 'MATERIAL_CITATION', MATERIAL_SAMPLE', 'PRESERVED_SPECIMEN', 'FOSSIL_SPECIMEN',
#' 'LIVING_SPECIMEN', 'LITERATURE', UNKNOWN' and 'OCCURRENCE'. Default setting removes
#' specimens and unknown observations.
#' Description may be found here: https://docs.gbif.org/course-data-use/en/basis-of-record.html, 
#' https://gbif.github.io/parsers/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html.
#' @param add_infos Character. Infos that may be added to the default output information.
#' List of IDs may be found at: https://techdocs.gbif.org/en/data-use/download-formats.
#' Default IDs contain 'taxonKey', 'scientificName', 'acceptedTaxonKey',
#' 'acceptedScientificName', 'individualCount', 'decimalLatitude', 'decimalLongitude',
#' 'basisOfRecord', 'coordinateUncertaintyInMeters', 'countryCode', 'country', 'year', 'datasetKey', 
#' 'institutionCode', 'publishingOrgKey', 'taxonomicStatus', 'taxonRank' and 'degreeOfEstablishment'. 
#' @param time_period Numerical vector. Observations will be downloaded according to the chosen
#' year range. Default is c(1000,3000). Observations with year = NA are kept by default.
#' @param identic_xy Logical. Should records with identical xy be kept?
#' @param wConverted_xy Logical. Should incorrectly converted lon/lat be kept?
#' Uses cd_ddmm() from 'CoordinateCleaner' R package.
#' @param centroids Logical. Should species records from raster centroids be kept?
#' Uses cd_round() from 'CoordinateCleaner' R package.
#' @param ntries Numeric. In case of failure from GBIF server or within the rgbif package, how many
#' download attempts should the function request? Default is '10' with a 2 seconds interval
#' between tries. If the attempts failed, an empty data.frame is return by default.
#' @param error_skip Logical. Should the search process continues if ntries failed ?
#' @param occ_samp Numeric. Determine how many GBIF occurrences will be sampled per geographic
#' tiles of the fragmented study area. Default is the maximum number of GBIF observations found
#' in a tile (i.e. ~100'000 records). A lower number may be set (<99'000) if the user only wants
#' a sample of the species GBIF observations, hence increasing the download process and the
#' generation of its range map if get_range() is employed afterwards.
#' @param ... Additonnal parameters for the function cd_round() of CoordinateCleaner.
#' @details Argument `grain` used for two distinct gbif records filtering. (1) Records filtering
#' according to gbif 'coordinateUncertaintyInMeters'; every records uncertainty > grain/2
#' are removed. Note: Records with no information on coordinate uncertainties are kept by
#' default. (2) Records filtering according to the number of longitude/latitude decimals;
#' if 110km < grain <= 11km, lon/lat with >= 1 decimal are kept, if 11km < grain <= 1100m,
#' lon/lat with >= 2 decimals kept; if 1100m < grain <= 110m, lon/lat with >= 3 decimals
#' are kept; if 110m < grain <= 11m, lon/lat with >= 4 decimals are kept;
#' if 11m < grain <= 1.1m, lon/lat with >= 5 decimals are kept etc...
#' @return Object of class data.frame with requested GBIF information. Although the function
#' works accurately, error outputs might still occur depending on the 'sp_name' used.
#' Therefore, default information detailed in 'add_infos' is stored so that sanity checks
#' may still be applied afterwards. Although crucial preliminary checks of species records
#' are done by the function, additional post exploration with the CoordinateCleaner R
#' package is still highly recommended.
#' @references 
#' Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann,
#' N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the
#' European Alps. Ecological monographs, 91(2), e01433. 10.1002/ecm.1433
#' 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity
#' information facility API. 10.5281/zenodo.6023735
#' 
#' Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli,
#' A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection
#' databases. Methods in Ecology and Evolution, 10(5), 744-751. 10.1111/2041-210X.13152
#' 
#' Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). Terra - CRAN
#' @seealso The (1) rgbif and (2) CoordinateCelaner packages for additional and more general
#' approaches on (1) downloading GBIF observations and (2) post-filtering those.
#' @example inst/examples/get_gbif_help.R
#' @importFrom terra ext
#' @importFrom rgbif name_backbone occ_data
#' @importFrom CoordinateCleaner cd_ddmm cd_round
#' @export
get_gbif = function(sp_name = NULL,
					search = TRUE,
					rank = NULL,
					phylum = NULL,
					class = NULL,
					order = NULL,
					family = NULL,
					conf_match = 80,
					geo = NULL,
					grain = 1000,
					duplicates = FALSE,
					absences = FALSE,
					no_xy = FALSE,
					basis =  c('OBSERVATION', 'HUMAN_OBSERVATION', 'MACHINE_OBSERVATION','OCCURRENCE',
						'MATERIAL_CITATION', 'MATERIAL_SAMPLE','LITERATURE'),
					add_infos = NULL,
					time_period = c(1000,3000),
					identic_xy = FALSE,
					wConverted_xy = FALSE,
					centroids = FALSE,
					ntries = 10,
					error_skip = TRUE,
					occ_samp = 99000,
					...) {


	######################################################
	### Stop message
	######################################################


	if (is.null(sp_name)|is.na(sp_name)) {
		stop("Given 'sp_name' is NA or NULL, a string must be provided...")
	}


	######################################################
	### Parameter config
	######################################################


	# For precision
	deci_preci = list(seq(0,10,1),rev(0.000011 * 10^(0:10)))
	deci_chosen = deci_preci[[1]][which(grain>=deci_preci[[2]])[1]]

	# For study area
	if (!is.null(geo)) {
		if (!(class(geo)%in%"SpatExtent")) {geo = terra::ext(geo)}
	} else {
		geo = terra::ext(-180,180,-90,90)
	}

	# For fields
	gbif_info = c('taxonKey','scientificName','acceptedTaxonKey','acceptedScientificName',
		'individualCount','occurrenceStatus','establishmentMeans','decimalLatitude','decimalLongitude',
		'basisOfRecord','coordinateUncertaintyInMeters','countryCode','country', 'year','datasetKey',
		'institutionCode','publishingOrgKey','taxonomicStatus','taxonRank',add_infos)
	gbif_info = gbif_info[order(gbif_info)]

		# For 'cd_ddmm'
	diff_records = list(c(0,5000,100000 * 10^(0:9)),rev(0.00000000001 * 10^(0:11)))
	xy_span = 10*10^(-(deci_chosen-1))

		# Create an empty ouptut
	e_output = data.frame(matrix(ncol=length(gbif_info),nrow=0,
		dimnames=list(NULL,gbif_info)))


	# Set geo.ref
	geo_ref = paste0("POLYGON((",geo$xmin," ",geo$ymin,", ",
							 	geo$xmax," ",geo$ymin,", ",
							  geo$xmax," ",geo$ymax,", ",
							  geo$xmin," ",geo$ymax,", ",
							  geo$xmin," ",geo$ymin,"))")


	######################################################
	### Backbone harmonization and N records
	######################################################


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

          if (any(bone_search$status%in%"ACCEPTED") & length(unique(key_test[,1]))==1){
            bone_search = bone_search[bone_search$status%in%"ACCEPTED",]
          }

        } else {
          bone_search = s_keep
        }
        # If not the same species overall return empty
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

	# Get the accepetedKey
	if (bone_search$status %in% "SYNONYM"){
		sp_key = bone_search$acceptedUsageKey
	} else {
		sp_key = bone_search$usageKey
	}

	# Check number of records in 'geo' first
	gbif_records = rgbif::occ_data(taxonKey=sp_key,limit=1,hasCoordinate=!no_xy,
			hasGeospatialIssue=FALSE,geometry=geo_ref)[1]$meta$count

	cat(">>>>>>>> Total number of records:",gbif_records,"\n")

	# Cancel request if n=0
	if (gbif_records==0) {
		cat("No species records found...","\n")
		return(e_output)
	}


	######################################################
	### Find the tiles where obs. < 100'000 observations
	######################################################


	## 1) If species records > 100'000, search for the optimum tiles
	if (gbif_records > 99000)
	{
		cat(">>>>>>>> Too many records: Retrieving relevant geographic tiles...","\n")

		# Start with 200 tiles
		tile_100 = make_tiles(geo,Ntiles=200,sext=TRUE)
		geo_tiles = tile_100[[1]]
		geo_meta = lapply(tile_100[[2]],function(x) x[])

		# Check number of records for each tiles
		gbif_tiles =
		sapply(geo_tiles, function(x) {
			gt_out = rgbif::occ_data(taxonKey=sp_key,limit=1,hasCoordinate=!no_xy,
				hasGeospatialIssue=FALSE,geometry=x)[1]$meta$count
			return(gt_out)
		})

		# We make sure that each tile will be fragmented enough
		# (i.e., < 100'000 species records each)
		keep_tiles =
		lapply(1:length(gbif_tiles),function(x){

			# Which tile
			tile_n = gbif_tiles[x]

			# Return NULL or POLYGON if 0 records or > 100'000 records
			if (tile_n==0) {return(NULL)}
			if (tile_n<99000) {return(geo_tiles[x])}

			# Keep fragmenting if > 100'000 records
			new_geo = list(geo_meta[[x]])
			pol_store=list()
			while (length(new_geo)!=0)
			{
				m_process =
				lapply(1:length(new_geo),function(y){

					# Convert to extent and create 100 new micro tiles
					new_ext = terra::ext(new_geo[[y]])
					micro_100 = make_tiles(new_ext,Ntiles=100,sext=TRUE)

					# If micro.100 is NULL return NULL (in case of an incorrect extent)
					if (is.null(micro_100)) {
						return(list(NULL,NULL))
					}

					# Continue
					micro_tiles = micro_100[[1]]
					micro_meta = lapply(micro_100[[2]],function(z) z[])

					# And count records again
					gbif_micro =
						sapply(micro_tiles, function(z) {
							gt_out = try(rgbif::occ_data(taxonKey=sp_key,limit=1,hasCoordinate=!no_xy,
							hasGeospatialIssue=FALSE,geometry=z)[1]$meta$count,silent=TRUE)
						return(gt_out)
					})

					# Transform in 0 when the polygon is too small for occ_data
					gbif_micro[sapply(gbif_micro,function(z) grepl("Error",z))]=0
					gbif_micro=as.numeric(gbif_micro)

					# Store if tiles with < 99'000 observations found
					goody = micro_tiles[gbif_micro<99000 & gbif_micro!=0]
					bady = micro_meta[!gbif_micro<99000 & gbif_micro!=0]
					return(list(goody,bady))
				})

				# Restructure
				all_good = lapply(m_process,function(y) y[[1]])
				all_good[sapply(all_good, is.null)] = NULL
				all_bad = lapply(m_process,function(y) y[[2]])
				all_bad[sapply(all_bad, is.null)] = NULL

				# Store if no blocks with > 99'000 observations found
				pol_store = c(pol_store,unlist(all_good,recursive=FALSE))
				new_geo = unlist(all_bad,recursive=FALSE)

				# Finally remove too small windows (for "invisible" GBIF duplicated records)
				inf_id = sapply(new_geo,function(y) y[2]-y[1]<1e-7)
				new_geo = new_geo[!inf_id]
			}
			return(unlist(pol_store,recursive=FALSE))
		})

		# Final list of POLYGONS with n observations < 100'000 observations
		geo_ref = unlist(keep_tiles)
	}
	
	######################################################
	#################### API Search ######################
	######################################################

	# Infromation messages
	if (occ_samp!=99000) {
		cat("...GBIF records of",sp_name,": download of sample of records starting...","\n")
	} else {
		cat("...GBIF records of",sp_name,": download of all records starting...","\n")
	}

	# Run the gbif search with the acceptedName per chosen tiles
	batch_search =
	lapply(1:length(geo_ref),function(x) {

		cat("\r",round(x*100/length(geo_ref),2),"%...")

		## Try the download first: may be request overload problems
		go_tile = geo_ref[x]
		gbif_search = try(
			rgbif::occ_data(taxonKey=sp_key,limit=occ_samp,hasCoordinate=!no_xy,
				hasGeospatialIssue=FALSE,geometry=go_tile),
		silent=TRUE)

		# If problems, just try to rerun with while with n attempts, otherwise return NULL
		if (class(gbif_search) %in% "try-error") {
			print(gbif_search[1])
			warning("GBIF query overload or rgbif package error [taxonKey=",sp_key,"]...","\n",sep="")

			# While
			if (class(gbif_search) %in% "try-error") {
				j = 0
				while (class(gbif_search) %in% "try-error" & j<ntries)
				{
					cat("Attempt",j+1,"...","\n")
					j = j+1
					gbif_search = try(
						rgbif::occ_data(taxonKey=sp_key,limit=occ_samp,hasCoordinate=!no_xy,
							hasGeospatialIssue=FALSE,geometry=go_tile),
					silent=TRUE)
					Sys.sleep(2)
				}
				if (class(gbif_search) %in% "try-error") {
					if (error_skip){
						cat("Attempts to download failed...Returning no results")
						return(e_output)
					} else {
						stop("ERROR (not skipped) for [taxonKey=",sp_key,"]...","\n",sep="")
					}
				}
			}
		}

		# If no results
		if (is.null(gbif_search$data)){
			
			return(e_output)
		
		} else {

			# Convert to a data.frame is needed
			if (class(gbif_search$data)[1]!="data.frame"){
				gbif_search = as.data.frame(gbif_search$data)
			}

			# Reordering data.frame and correcting if missing columns
			gbif_reorder = gbif_search[,order(names(gbif_search))]
			missing_col = gbif_info[!(gbif_info%in%names(gbif_reorder))]
			gbif_reorder[,missing_col] = NA
			gbif_correct = gbif_reorder[,order(names(gbif_reorder))]

			# Remove row we don't want
			return(gbif_correct[,gbif_info])
		}
	})

	# Combine all results in one data.frame
	gbif_compile = do.call("rbind",batch_search)

	# Keep specific fields (just for safety)
	gbif_compile = gbif_compile[,gbif_info]


	#############################################################
	#################### Records filtering ######################
	#############################################################


	############ 1) Grain filtering
	cat("\n","---> Grain filtering...","\n",sep="")
		
		# GBIF uncertainty
	id_certain = gbif_compile$coordinateUncertaintyInMeters<=grain/2
	id_certain[is.na(id_certain)] = TRUE
	gbif_correct = gbif_compile[id_certain,]

		# GBIF lon/lat decimals
	if (grain < 1.1e5) {

		# Remove latitude or/and longitude with no decimals
		lonlat_format = data.frame(decimalLatitude=as.character(gbif_correct$decimalLatitude),
			decimalLongitude=as.character(gbif_correct$decimalLongitude))
		id_deci = grepl("\\.",lonlat_format[,1]) + grepl("\\.",lonlat_format[,2])
		lonlat_deci = lonlat_format[id_deci %in% 2,]
		gbif_correct = gbif_correct[id_deci %in% 2,]

		# Keep coordinates compatible with the input 'grain'
		declat = gsub(".*\\.","",lonlat_deci[,1])
		declon = gsub(".*\\.","",lonlat_deci[,2])
		id_grain = nchar(declon)>=deci_chosen & nchar(declat)>=deci_chosen
		gbif_correct = gbif_correct[id_grain,]

	} else {
		id_grain = NULL
	}

		# Removal summary
	if (any(names(table(c(id_certain,id_grain))) %in% FALSE)){
		removed = table(c(id_certain,id_grain))[1]
		cat("Records removed:",removed,"\n")
	} else {
		cat("Records removed:",0,"\n")
	}


	############ 2) Removing xy duplicates
	if (!duplicates){

		cat("---> Removal of duplicated records...","\n")
		
		id_dup = !duplicated(gbif_correct[,c("decimalLongitude","decimalLatitude")])
		gbif_correct = gbif_correct[id_dup,]

		# Removal summary
		if (any(names(table(id_dup))%in%FALSE)){
			removed = table(id_dup)[1]
			cat("Records removed:",removed,"\n")
		} else {
			cat("Records removed:",0,"\n")
		}
	}

	
	############ 3) Removing absences
	if (!absences){

		cat("---> Removal of absence records...","\n")
		
		id_abs = !(gbif_correct$individualCount %in% 0 | gbif_correct$occurrenceStatus %in% "ABSENT")
		gbif_correct = gbif_correct[id_abs,]

		# Removal summary
		if (any(names(table(id_abs))%in%FALSE)){
			removed = table(id_abs)[1]
			cat("Records removed:",removed,"\n")
		} else {
			cat("Records removed:",0,"\n")
		}
	}


	############ 4) Select basis of records
	cat("---> Basis of records selection...","\n")

	id_basis = gbif_correct$basisOfRecord %in% basis
	gbif_correct = gbif_correct[id_basis,]

		# Removal summary
	if (any(names(table(id_basis))%in%FALSE)){
		removed = table(id_basis)[1]
		cat("Records removed:",removed,"\n")
	} else {
		cat("Records removed:",0,"\n")
	}

	
	############ 5) Select records according to year range
	cat("---> Time period selection...","\n")

	id_year = gbif_correct$year >= min(time_period) & gbif_correct$year <= max(time_period)
	id_year[is.na(id_year)] = TRUE
	gbif_correct = gbif_correct[id_year,]

		# Removal summary
	if (any(names(table(id_year))%in%FALSE)){
		removed = table(id_year)[1]
		cat("Records removed:",removed,"\n")
	} else {
		cat("Records removed:",0,"\n")
	}

	
	############ 6) Remove records with identical xy
	if (nrow(gbif_correct) > 0){
		if (!identic_xy){
		
			cat("---> Removal of identical xy records...","\n")

			id_diff = c(!(abs(gbif_correct[,"decimalLatitude"]) ==
				abs(gbif_correct[,"decimalLongitude"])))
			gbif_correct = gbif_correct[id_diff,]

			# Removal summary
			if (any(names(table(id_diff))%in%FALSE)){
				removed = table(id_diff)[1]
				cat("Records removed:",removed,"\n")
			} else {
				cat("Records removed:",0,"\n")
			}
		}
	}


	############ 7) Remove wrongly lon/lat converted xy (only if > 50 records in one dataset)
	gbif_correct = as.data.frame(gbif_correct)
	gbif_nrow = nrow(gbif_correct)

	if (nrow(gbif_correct) > 0){
		if (!wConverted_xy){

			cat("---> Removal of wrong lon/lat converted records...","\n")

			# Keep dataset with NAs
			gbif_na = gbif_correct[is.na(gbif_correct$datasetKey),]

			# Summary of datasets & cd_ddmm parameter choice
			gbif_datasets = names(table(gbif_correct$datasetKey))
			if (nrow(gbif_correct) < 10000) {mat_size = 100} else {mat_size = 1000}

			# Apply correction if the dataset > 50 records
			gbif_ddmm =
			lapply(1:length(gbif_datasets), function(x){

				gbif_dataset = gbif_correct[gbif_correct$datasetKey%in%gbif_datasets[x],]
				if (nrow(gbif_dataset) > 50){
					diff_chosen = diff_records[[2]][which(nrow(gbif_dataset)<diff_records[[1]])[1]]
					gbif_dataset = CoordinateCleaner::cd_ddmm(gbif_dataset,lon="decimalLongitude",
						lat="decimalLatitude",ds="datasetKey",mat_size=mat_size,diff=diff_chosen,min_span=xy_span)
				}
				return(gbif_dataset)
			})
			gbif_Wna = do.call("rbind", gbif_ddmm)
			gbif_correct = rbind(gbif_na,gbif_Wna)

			# Removal summary
			cat("Records removed:",gbif_nrow-nrow(gbif_correct),"\n")
		}
	}


	############ 8) Remove raster centroids (only if > 100 records in one dataset) (ref:coordinate cleaner article in methods)
	gbif_nrow = nrow(gbif_correct)

	if (nrow(gbif_correct) > 0){
		if (!centroids){

			cat("---> Removal of raster centroids...","\n")

			# Keep dataset with NAs
			gbif_na = gbif_correct[is.na(gbif_correct$datasetKey),]
			
			# Summary of datasets & cd_ddmm parameter choice
			gbif_datasets = names(table(gbif_correct$datasetKey))

			# Apply correction if the dataset > 100 records
			gbif_round =
			lapply(1:length(gbif_datasets), function(x){

				gbif_dataset = gbif_correct[gbif_correct$datasetKey%in%gbif_datasets[x],]
				if (nrow(gbif_dataset) > 100){
					gbif_temp = try(CoordinateCleaner::cd_round(gbif_dataset,lon="decimalLongitude",
						lat="decimalLatitude",ds="datasetKey",graphs=FALSE,...),silent=TRUE)
					if (class(gbif_temp)%in%"try-error") {
						return(gbif_dataset)
					} else {
						gbif_dataset=gbif_temp
					}
				}
				return(gbif_dataset)
			})
			gbif_Wna = do.call("rbind", gbif_round)
			gbif_correct = rbind(gbif_na,gbif_Wna)

			# Removal summary
			cat("Records removed:",gbif_nrow-nrow(gbif_correct),"\n")
		}
	}
	if (nrow(gbif_correct)==0) {
		cat("No records left after selection...","\n")
		return(e_output)
	} else {
		return(cbind(input.search=sp_name,gbif_correct))
	}
}
