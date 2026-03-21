### ==================================================================
### evaluate_range
### ==================================================================
#' Evaluate Range Maps Against Independent Validation Data
#' 
#' Compare range maps produced by \code{get_range()} with external validation
#' data, such as species distribution models (SDMs) or expert-derived range
#' maps, and summarize precision, sensitivity, specificity, and TSS.
#'
#' @param root_dir Character string giving the root directory that contains both
#' the generated range maps and the validation data.
#' @param valData_dir Character string giving the directory, relative to
#' \code{root_dir}, containing the validation data.
#' @param ecoRM_dir Character string giving the directory, relative to
#' \code{root_dir}, containing range maps generated with \code{get_range()}.
#' @param valData_type Character string indicating the expected validation-data
#' format: \code{"SHP"} or \code{"TIFF"}.
#' @param verbose Logical. Should progress information be printed while the
#' comparison is running?
#' @param print_map Logical. If \code{TRUE}, write a PDF overlay map for each
#' evaluated species.
#' @param mask Optional \code{SpatRaster} used as a study-area mask and common
#' comparison domain.
#' @param res_fact Optional integer aggregation factor used to coarsen the
#' native resolution before comparison.
#' @details TIFF validation files must have file names matching the species
#' names of the range maps. Shapefile-based validation data must include a
#' column named \code{sci_name} with matching species names.
#'
#' The function can optionally mask the focal study region and aggregate maps to
#' coarser resolutions before calculating evaluation metrics, which is useful
#' when comparing products with different native resolutions.
#' @return A list with two elements: \code{df_eval}, a data frame containing
#' per-species evaluation statistics, and \code{overlay_list}, a list of raster
#' overlays used for plotting and inspection.
#' @references
#' Pinkert, S., Sica, Y. V., Winner, K., & Jetz, W. (2023). The potential of 
#' ecoregional range maps for boosting taxonomic coverage in ecology and 
#' conservation. Ecography, 12, e06794.
#' @example inst/examples/evaluate_range_help.R
#' @importFrom terra rast ext crs project aggregate rasterize crop extend
#' resample values classify ncol nrow plot
#' @importFrom sf st_read st_as_sf st_union st_drop_geometry st_transform
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics legend par
#' @export
evaluate_range <- function(root_dir = NULL,
                       valData_dir = NULL,
                       ecoRM_dir = NULL,
                       valData_type = NULL,
                       verbose = TRUE,
                       print_map = TRUE,
                       mask = NULL,
                       res_fact = NULL) {

  ######################################################
  ### Stop messages
  ######################################################


  check_null_na(root_dir, "root_dir")
  check_null_na(valData_dir, "valData_dir")
  check_null_na(ecoRM_dir, "ecoRM_dir")
  if (!valData_type%in%c("SHP","TIFF")) {
    stop("Given 'valData_type' must be SHP or TIFF...")
  }
  check_null_na(verbose, "verbose")
  check_null_na(print_map, "print_map")


  ######################################################
  ### Code
  #####################################################


  # No root_dir
  if (is.null(root_dir)){
    stop(
      "Please provide a main working directory 'root_dir'... \n"
      )
  }

  # Create dummy objects 
  overlay.list <- list()
  df.eval <- NULL
  
  # Get ecoRM file list
  f.list.eco.full <- list.files(
    path = file.path(root_dir, ecoRM_dir),
    pattern = ".tif",
    full.names = TRUE
  )
  f.list.eco <- sub("\\.tif$", "", basename(f.list.eco.full))
  
  # Retrieve valData
  shp.names <- unique(
    list.files(
      path = file.path(root_dir, valData_dir),
      pattern = "\\.shp$",
      full.names = TRUE
    )
  )
  
  if (length(shp.names) != 0) {
    shp.files <- do.call(rbind, lapply(shp.names, sf::st_read))

    if (any(duplicated(shp.files$sci_name))) {
      unique.sci.names <- unique(shp.files$sci_name)

      grouped.geometries <-
      lapply(unique.sci.names, function(name) {
        geoms <- shp.files[sf::st_drop_geometry(shp.files)$sci_name == name, ]
        sf::st_union(geoms)
      })

      shp.files <- data.frame(
        sci_name = unique.sci.names,
        geometry = do.call(c, grouped.geometries)
      )
      shp.files <- sf::st_as_sf(shp.files)
    }

    f.list.valRM <- shp.files$sci_name

  } else {
    f.list.valRM.full <- list.files(
      path = file.path(root_dir, valData_dir),
      pattern = "\\.tif$",
      full.names = TRUE,
      recursive = TRUE
    )

    f.list.valRM <- sub("\\.tif$", "", basename(f.list.valRM.full))
    valData_type <- "TIFF"
  }
  
  # Check and report match of file names
  if (length(f.list.valRM) == 0) {
    stop(paste(
      "No names of 'valData_dir' match those of the ecoRMs.",
      "Assuming that your species names in the ecoRMs are separated",
      "by spaces and followed by '.tif'. Text before the species",
      "names (all separated with underscores) is ignored."
    ))
  }
  
  cat("-Note-", "\n")
  if (valData_type == "TIFF") {
    f.list.matches <- intersect(f.list.eco, basename(f.list.valRM))
    cat(
      sprintf(
        paste(
          "%.2f%% (%d) of the species names of ecoregions",
          "match with names of the validation data files"
        ),
        100 * length(f.list.matches) / length(f.list.valRM),
        length(f.list.matches)
      )
    )

  } else {
    f.list.matches <- intersect(f.list.eco, f.list.valRM)
    cat(
      sprintf(
        paste(
          "%.2f%% (%d) of the species names of ecoregions",
          "match with those in the validation data column 'sci_name'"
        ),
        100 * length(f.list.matches) / length(f.list.valRM),
        length(f.list.matches)
      )
    )
  }

  cat("------", "\n")
  
  df.eval <- data.frame(
    species = f.list.matches,
    nbr_pres = NA,
    nbr_pres_ecoRM = NA,
    nbr_true_pres = NA,
    nbr_false_pres = NA,
    Prec_ecoRM = NA,
    Sen_ecoRM = NA,
    Spec_ecoRM = NA,
    TSS_ecoRM = NA,
    type = ecoRM_dir
  )
  
  # Process each species
  for (i in seq_along(f.list.matches)) {
    if (verbose) cat(i, " Species: ", f.list.matches[i], "\n")
    
    if (is.null(mask)) {
      domain.raster <- NULL
    } else {
      domain.raster <- mask
    }
    
    # Load ecoRM raster
    ecoRM <- terra::rast(
      f.list.eco.full[grep(f.list.matches[i], f.list.eco.full)]
    )
    
    # Load validation data
    if (valData_type == "TIFF") {
      valRM <- terra::rast(
        f.list.valRM.full[grep(f.list.matches[i], f.list.valRM.full)]
      )
      suppressWarnings({
        if (isTRUE(all.equal(terra::crs(ecoRM), terra::crs(valRM)))==FALSE) {
          valRM <- terra::project(valRM, ecoRM)
        }
      })

    } else {
      valRM <- shp.files[grep(f.list.matches[i], shp.files$sci_name), ]
      valRM <- sf::st_transform(valRM, crs = terra::crs(ecoRM))
    }
    
    ext.e <- terra::ext(ecoRM)
    ext.v <- terra::ext(valRM)
    
    if (is.null(mask)) {

      if (isTRUE(all.equal(ext.e, ext.v))==FALSE) {
        
        combined.extent <- terra::ext(
          min(ext.e[1], ext.v[1]),
          max(ext.e[2], ext.v[2]),
          min(ext.e[3], ext.v[3]),
          max(ext.e[4], ext.v[4])
        )
        domain.raster <- terra::extend(ecoRM, combined.extent)

        if (valData_type == "SHP") {
          valRM <- terra::rasterize(
            x = valRM,
            y = domain.raster,
            field = 1,
            background = NA
          )
        }

      } else {
        if (valData_type == "SHP") {
          valRM <- terra::rasterize(
            x = valRM,
            y = ecoRM,
            field = 1,
            background = NA
          )
        }
        domain.raster <- valRM
      }

      domain.raster[] <- 0
      domain.raster[which(is.na(valRM[]))] <- NA
    }
    
    # Resample domain raster if resolution factor is provided
    if (!is.null(res_fact) && res_fact >= 1) {
      domain.raster <- terra::aggregate(
        x = domain.raster,
        fact = res_fact,
        fun = max,
        na.rm = TRUE
      )
    }
    
    # Ensure resolution and extent match
    ecoRM <- terra::resample(ecoRM, domain.raster, method = "near")
    valRM <- terra::resample(valRM, domain.raster, method = "near")
    
    if (!is.null(mask)) {
      ecoRM[is.na(domain.raster[])] <- NA
      valRM[is.na(domain.raster[])] <- NA
    }
    
    # Create a binarized ecoRM
    ecoRM.bin <- ecoRM
    ecoRM.bin[!is.na(domain.raster[])] <- 0
    ecoRM.bin[ecoRM == 1] <- 1
    
    # Reclassify validation raster
    if (all(range(terra::values(valRM), na.rm = TRUE) == c(0, 1))) {
      valRM.pres <- terra::classify(
        x = valRM,
        rcl = rbind(c(0, 0.5, 0),c(0.5, 1, 1)),
        include.lowest = TRUE
      )
    
    } else {
      valRM.pres <- valRM
      valRM.pres[!is.na(domain.raster[])] <- 0
      valRM.pres[valRM == 1] <- 1
    }
    
    # Burn into the overlay layer with domain_raster properties
    overlay.raster <- valRM.pres * 2 + ecoRM.bin
    overlay.raster[is.na(domain.raster[])] <- NA 
    
    names(overlay.raster) <- f.list.matches[i]
    overlay.list[[i]] <- overlay.raster
    
    # Evaluate presence and absence
    df.eval$nbr_pres[i] <- sum(overlay.raster[] %in% c(2,3), na.rm = TRUE)
    df.eval$nbr_pres_ecoRM[i] <- sum(overlay.raster[] %in% c(1,3), na.rm = TRUE)
    df.eval$nbr_true_pres[i] <- sum(overlay.raster[] == 3, na.rm = TRUE)
    df.eval$nbr_false_pres[i] <- sum(overlay.raster[] == 1, na.rm = TRUE)
    nbr.abs <- sum(overlay.raster[] %in% c(0, 1), na.rm = TRUE)
    nbr.true.abs <- sum(overlay.raster[] == 0, na.rm = TRUE)

    df.eval$Prec_ecoRM[i] <-
      df.eval$nbr_true_pres[i] /
        (df.eval$nbr_true_pres[i] + df.eval$nbr_false_pres[i])
    df.eval$Sen_ecoRM[i] <-
      df.eval$nbr_true_pres[i] / df.eval$nbr_pres[i]
    df.eval$Spec_ecoRM[i] <-
      nbr.true.abs / nbr.abs
    df.eval$TSS_ecoRM[i] <-
      df.eval$Sen_ecoRM[i] + df.eval$Spec_ecoRM[i] - 1
    
    # Plot the overlay raster
    if (verbose && print_map) {
      
      if (!dir.exists(file.path(root_dir, "eval_output"))) {
        dir.create(file.path(root_dir, "eval_output"))
      }
      
      aspect.r <- terra::nrow(overlay.raster) / terra::ncol(overlay.raster)
      colors <- c("gray", "red", "blue", "purple")
      breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
      
      pdf(file = file.path(root_dir,"eval_output",
        paste0("Evaluation_map_", f.list.matches[i], ".pdf")),
        height = (11 * aspect.r)+3, width = 11
      )
      par(mfrow = c(1, 1))
      terra::plot(
        overlay.raster,
        col = colors,
        breaks = breaks,
        legend = FALSE,
        main = paste0(
          "Overlay validation RM and ecoRM - ", "Species: ",
          f.list.matches[i],"\n", 
          "Precision (TP/[TP+FP]) = ",
          round(df.eval$Prec_ecoRM[i], digits = 2)," & ", 
          "Sensitivity (TP/[TP+FA]) = ",
          round(df.eval$Sen_ecoRM[i], digits = 2)),
        las = 1
      )
      
      # Adding a legend
      legend("bottomright",
        legend = c(
          "Abs in both (TA)",
          "Pres in ecoRM only (FP)",
          "Pres in valRM only (FA)",
          "Pres in both (TP)"
        ),
        fill = colors,
        bg = NA,
        box.col = NA,
        inset = c(0,0.1)
      )
      dev.off()
    }
  }
  
  if (verbose) {
    cat("Cross-species mean Prec & Sensitivity:",
      round(
        mean(rowMeans(df.eval[, c("Prec_ecoRM", "Sen_ecoRM")],na.rm = TRUE)),
        digits = 2),
      "\n"
    )
  }
  
  if (print_map) {
    cat("### Maps have been saved to:",
      file.path(root_dir, "eval_output"), "###\n")
  }
  
  output <- list(df_eval = df.eval, overlay_list = overlay.list)
  return(output)
}
