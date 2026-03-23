### ==================================================================
### get_status
### ==================================================================
#' Retrieve Taxonomy and IUCN Status from GBIF
#'
#' Query the GBIF backbone taxonomy for a species name, report the matched
#' accepted name and synonyms, and return the associated IUCN status when
#' available.
#'
#' Use \code{get_status()} before \code{get_gbif()} when you want to inspect
#' how GBIF resolves an input name, which name is treated as the accepted
#' backbone taxon, and which synonyms are included in the taxon concept used
#' for occurrence retrieval.
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
#' @param all Logical. If \code{FALSE} (default), return the accepted name and
#' its synonyms. If \code{TRUE}, also include children and related names that
#' are not used by \code{get_gbif()}.
#' @details When \code{all = FALSE}, the returned rows correspond to the
#' accepted name and synonyms linked to the accepted GBIF taxon key. This is
#' the same taxon concept used internally by \code{get_gbif()}.
#'
#' When \code{all = TRUE}, the output is expanded with \code{CHILDREN} and
#' \code{RELATED} names returned by the GBIF backbone. These extra rows are
#' useful for taxonomic inspection, but they are not included when
#' \code{get_gbif()} downloads occurrence records.
#' @return A data frame with the columns \code{canonicalName}, \code{rank},
#' \code{gbif_key}, \code{scientificName}, \code{gbif_status},
#' \code{Genus}, \code{Family}, \code{Order}, \code{Class},
#' \code{Phylum}, \code{IUCN_status}, and \code{sp_nameMatch}. The row with
#' \code{gbif_status = "ACCEPTED"} identifies the accepted GBIF taxon concept
#' used by \code{get_gbif()}, while \code{sp_nameMatch} marks the row that most
#' closely matches the submitted input name.
#' @references 
#' Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the
#' global biodiversity information facility API. 10.5281/zenodo.6023735
#' @seealso \code{get_gbif()} to download occurrences for the accepted taxon
#' concept returned here; the \code{rgbif} package for more general approaches
#' to querying the GBIF backbone taxonomy.
#' @example inst/examples/get_status_help.R
#' @importFrom rgbif name_backbone name_usage
#' @importFrom methods is
#' @export
get_status <- function(sp_name = NULL,
                    search = TRUE,
                    rank = NULL,
                    phylum = NULL,
                    class = NULL,
                    order = NULL,
                    family = NULL,
                    conf_match = 80,
                    all = FALSE)
{

  ######################################################
  ### Stop messages
  ######################################################


  # General
  check_character_vector(sp_name, "sp_name")
  check_logical(search, "search")
  check_numeric(conf_match, "conf_match")
  check_logical(all, "all")


  ######################################################
  ### Get the status
  ######################################################


  # Empty output
  e.output <- data.frame(canonicalName = NA,
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
    bone.search <- rgbif::name_backbone(sp_name,
                                     rank = rank,
                                     phylum = phylum,
                                     class = class,
                                     order = order,
                                     family = family,
                                     verbose = FALSE,
                                     strict = FALSE)
  } else {
    # Search input name via strict match and refined search
    bone.search <- rgbif::name_backbone(sp_name,
                                       verbose = TRUE,
                                       strict = TRUE)

    q.crit <- !vapply(
      list(rank, phylum, class, order, family),
      is.null,
      FUN.VALUE = logical(1)
    )

    # Filter by the supplied criteria if results are available
    if (!bone.search$matchType[1] %in% "NONE"){ 

      if (any(q.crit)){
        id.crit <- c("rank", "phylum", "class", "order", "family")[q.crit]
        p.crit <- unlist(list(rank, phylum, class, order, family)[q.crit])
        n.test <- id.crit %in% names(bone.search)

        if (any(n.test)){
          # selecting which
          id.crit2 <- id.crit[n.test]
          p.crit2 <- p.crit[n.test]

          # Apply the requested criteria
          for (i in seq_along(id.crit2)){
            bone.search <-
              bone.search[c(bone.search[, id.crit2[i]])[[1]] %in% p.crit2[i], ]

            if (nrow(bone.search) == 0){
              bone.search <- data.frame(matchType = "NONE")
            }
          }
        }

        if (!base::all(n.test)){
          pp <- paste(id.crit[!n.test],collapse = ", ")
          warning(
            paste0(
              "'",pp,"' level(s) not available for this taxa in GBIF",
              " could not be employed..."
            )
          )
        }
      }
    }
    
    # Continue with or without the supplied criteria
    if (nrow(bone.search) > 1){

      if (base::all(!bone.search$rank %in% c("SPECIES", "SUBSPECIES", "VARIETY"))){
        cat("Not match found...", "\n")
        return(e.output)

      } else {
        s.keep <- bone.search[bone.search$rank
                      %in% c("SPECIES", "SUBSPECIES", "VARIETY"),]
        s.keep <- s.keep[s.keep$status %in% c("ACCEPTED", "SYNONYM"),]
        if (nrow(s.keep) == 0){
          cat("Not match found...","\n")
          return(e.output)

        } else if (nrow(s.keep) > 1){

          # If only subspecies and variety are found,
          # we need to (default) prioritize
          if (base::all(s.keep$rank %in% c("VARIETY", "SUBSPECIES"))){

            if ("var." %in% strsplit(sp_name," ")[[1]]){
              bone.search <- s.keep[s.keep$rank %in% "VARIETY",]

              if (nrow(bone.search) == 0){
                bone.search <- s.keep[s.keep$rank %in% "SUBSPECIES",]
              }

            } else {
              bone.search <- s.keep[s.keep$rank %in% "SUBSPECIES",]
            }

          } else {
            bone.search <- s.keep[s.keep$rank %in% "SPECIES",]
          }

          focp <- c("familyKey", "orderKey", "classKey", "phylumKey")
          coltax <- focp %in% colnames(bone.search)
          key.test <- bone.search[, focp[coltax]]

          if (any(bone.search$status %in% "ACCEPTED")
                & length(unique(key.test[, 1][[1]])) == 1) {
            bone.search <- bone.search[bone.search$status %in% "ACCEPTED", ]
          }

        } else {
          bone.search <- s.keep
        }

        # If not the same species overall return NULL
        s.usp <- length(unique(bone.search$speciesKey)) == 1

        if (!s.usp){
          cat(
            paste(
              "No synonyms distinction could be made.",
              "Consider using phylum/class/order/family..."
            ),
          "\n")
          return(e.output)

        } else {
          bone.search <- bone.search[1, ]
        } 
      }
    }
  }

  if (bone.search$matchType %in% "NONE") {
    cat("No species name found...", "\n")
    return(e.output)
  }

  if (bone.search$confidence[1] < conf_match) {
    cat("Confidence match not high enough...", "\n")
    return(e.output)
  }  

   # Keep original scientific name
  sc.name <- suppressWarnings(bone.search$scientificName)

  # Extract key of accepted name
  if (bone.search$status %in% "SYNONYM") {
    accep.key <- bone.search$acceptedUsageKey

  } else {
    accep.key <- bone.search$usageKey
  }

  # Extract accepted name and save it with its key in the prepared output
  accep.name <- rgbif::name_usage(accep.key, data = "name")$data
  syn.syn <- rgbif::name_usage(accep.key, data = "synonyms")$data
  main.dat <-  rgbif::name_usage(accep.key, data = "all")$data
  iucn <- try(
    rgbif::name_usage(
      accep.key, data = "iucnRedListCategory"
    )$data,
    silent = TRUE
  )

  if (methods::is(iucn,"try-error")) {
    iucn <- "NOT_FOUND"

    } else {
      iucn <- iucn$category
    }

  # Avoid the canonicalName error (sometimes the column is absent wtf...)
  if (!"canonicalName" %in% colnames(main.dat)){
    accep.name$canonicalName <- NA
    syn.syn$canonicalName <- NA
    main.dat$canonicalName <- NA
  }

  # Specific columns
  c.key <- suppressWarnings(c(accep.key, syn.syn$key))
  c.sc <- suppressWarnings(
    c(accep.name$scientificName,
      syn.syn$scientificName)
  )
  c.can <- suppressWarnings(
    c(accep.name$canonicalName,
      syn.syn$canonicalName)
  )
  c.status <- c("ACCEPTED",
    rep("SYNONYM",
      length(suppressWarnings(syn.syn$scientificName)))
  )

  # If all=TRUE, then continue the search to find possible name correspondence
  if (all) {
    # Combine everything and search for related names
    all.key <- suppressWarnings(c(accep.key, syn.syn$key, main.dat$key))
    all.version <-
    lapply(all.key, function(x){
      
      out <- suppressWarnings(rgbif::name_usage(x, data = "related")$data)

      if (is.null(out)|nrow(out) == 0){
        return(NULL)
      
      } else {
        n.ref <- c("canonicalName", "scientificName")
        r.col <- n.ref[n.ref %in% names(out)]
        r.out <- data.frame(
          canonicalName = rep(NA,nrow(out)),
          scientificName = rep(NA,nrow(out))
        )
        r.out[,r.col] <- out[, r.col]
        return(
          data.frame(
            key = x,
            canonicalName = r.out$canonicalName,
            scientificName = r.out$scientificName
          )
        )
      }
    })

    # Extract all names
    accep.n <- suppressWarnings(
      accep.name[, c("canonicalName", "key", "scientificName")]
    )
    accep.n$key <- accep.key
    c.n <- suppressWarnings(
      main.dat[, c("canonicalName", "key", "scientificName")]
    )
    r.n <- suppressWarnings(
      unique(do.call("rbind", all.version))
    )

    # Conditions for synonymy
    syn.n <- try(
      suppressWarnings(
        syn.syn[, c("canonicalName", "key", "scientificName")]
      ), 
      silent = TRUE
    )
    if (class(syn.n)[1] %in% "try-error") {
      syn.n <- data.frame(
        canonicalName = NULL,
        key = NULL,
        scientificName = NULL
      )
    }
    all.names <- rbind(accep.n, syn.n, c.n, r.n)

    # Specific columns
    c.key <- all.names$key
    c.sc <- suppressWarnings(all.names$scientificName)
    c.can <- suppressWarnings(all.names$canonicalName)
    c.status <- c("ACCEPTED",
              rep("SYNONYM", nrow(syn.n)),
              rep("CHILDREN", nrow(c.n)),
              rep("RELATED", nrow(r.n)))
  }

  # Which is null in main.dat for Genus, Family, Order, Phyllum?
  gfocp <- c("genus", "family", "order", "class", "phylum")
  exist.not <- gfocp %in% names(main.dat)
  main.out <- data.frame(
    Genus = NA,
    Family = NA,
    Order = NA,
    Class = NA,
    Phylum = NA
  )
  main.out[exist.not] <- main.dat[, gfocp[exist.not]]

  # Extract accepted names and synonyms
  e.output <- data.frame(
    canonicalName = c.can,
    rank = bone.search$rank,
    gbif_key = c.key,
    scientificName = c.sc,
    gbif_status = c.status,
    main.out,
    IUCN_status = iucn,
    sp_nameMatch = bone.search$matchType
  )

  # Remove duplicated & # change sp_nameMatch to know which
  # row related to the sp_name input
  e.output <- e.output[!duplicated(e.output$scientificName), ]
  e.output[e.output$scientificName %in% sc.name, "sp_nameMatch"][1] <- "INPUT"

  # In case missing match
  if (!"INPUT"%in%e.output$sp_nameMatch) {
    
    # Normalize hybrid symbol spacing
    # (replace plain x or × with no spaces around)
    e.output$scientificName <-
            gsub("\\s*[x\u00D7]\\s*", "\u00D7", e.output$scientificName)
    sc.name <-
            gsub("\\s*[x\u00D7]\\s*", "\u00D7", sc.name)

    # Normalize spaces and trim
    e.output$scientificName <-
            gsub("\\s+", " ", trimws(e.output$scientificName))
    sc.name <-
            gsub("\\s+", " ", trimws(sc.name))

    # Use agrep for approximate matching
    # with max.distance of 0.1 (1% difference)
    matches <- agrep(
      sc.name,
      e.output$scientificName,
      max.distance = 0.01,
      ignore.case = TRUE
    )

    # Assign "INPUT" to the first matched index if any found
    if (length(matches) > 0) {
      e.output$sp_nameMatch[matches[1]] <- "INPUT"
    }
  }
  return(e.output)
}
