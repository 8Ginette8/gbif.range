### =========================================================================
### Function: evaluate the sensitivity & precision of range maps
### =========================================================================
#' Evaluates the sensitivity & precision of range maps based on validation
#' data, such as predictions of species distributions (SDMs) or IUCN expert
#' range maps. See Pinkert et al. (2023) for metrics and comparisions of 
#' various types of range data, including expert range maps, SDMs and ecoregional
#' range maps. Optional functionalities include the masking of the focal study 
#' region (see 'mask') and aggregations of the input maps to different resolutions,
#' given the importance of these factors for specific applications (Pinkert et al., 2023).
#' @param root.dir Character. Root directory to files
#' @param valData.dir Numeric. Buffer width parameter
#' @param ecoRM.dir  Numeric. Number of observation points
#' @param valData.type Character. Type of valData - either "SHP" or "TIFF"
#' @param verbose Logical. Optional - report details while running
#' @param print.map Logical. Optional - if verbose=TRUE should a overlap map be printed
#' @param mask ?
#' @param res.fact Integer. Factor for coarsening the original resolution
#' @references
#' Pinkert, S., Sica, Y. V., Winner, K., & Jetz, W. (2023). The potential of 
#' ecoregional range maps for boosting taxonomic coverage in ecology and 
#' conservation. Ecography, 12, e06794.
#' @export
#' @importFrom terra rast ext crs project aggregate rasterize crop extend resample values classify ncol nrow plot
#' @importFrom sf st_read st_as_sf st_union st_drop_geometry st_transform
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics legend par
#' @example inst/examples/evaluateRM_help.R

evaluateRM <- function(root.dir = NULL,
                       valData.dir = NULL,
                       ecoRM.dir = NULL,
                       valData.type = NULL,
                       verbose = TRUE,
                       print.map = TRUE,
                       mask = NULL,
                       res.fact = NULL) {
  # Create dummy objects 
  overlay_list <- list()
  df_eval <- NULL
  
  # Get ecoRM file list
  f_list_eco_full <- list.files(file.path(root.dir, ecoRM.dir),
                                pattern = ".tif",
                                full.names = TRUE)
  f_list_eco <- sub("\\.tif$", "", basename(f_list_eco_full))
  
  # Retrieve valData
  shp_names <- unique(list.files(
    file.path(root.dir, valData.dir),
    pattern = "\\.shp$",
    full.names = TRUE
  ))
  
  if (length(shp_names) != 0) {
    shp_files <- do.call(rbind, lapply(shp_names, sf::st_read))
    if (any(duplicated(shp_files$sci_name))) {
      unique_sci_names <- unique(shp_files$sci_name)
      grouped_geometries <- lapply(unique_sci_names, function(name) {
        geoms <- shp_files[sf::st_drop_geometry(shp_files)$sci_name == name, ]
        sf::st_union(geoms)
      })
      shp_files <- data.frame(sci_name = unique_sci_names, geometry = do.call(c, grouped_geometries))
      shp_files <- sf::st_as_sf(shp_files)
    }
    f.list_valRM <- shp_files$sci_name
  } else {
    f.list_valRM_full <- list.files(
      file.path(root.dir, valData.dir),
      pattern = "\\.tif$",
      full.names = TRUE,
      recursive = TRUE
    )
    f.list_valRM <- sub("\\.tif$", "", basename(f.list_valRM_full))
    valData.type = "TIFF"
  }
  
  # Check and report match of file names
  if (length(f.list_valRM) == 0) {
    stop(
      "No names of 'valData.dir' match those of the ecoRMs.
         Assuming that your species names of the ecoRMs are separated by a space and followed by '.tif'.
         Text before the species names (all separated with underscores) are ignored. \n"
    )
  }
  
  cat("-Note-", "\n")
  if (valData.type == "TIFF") {
    f.list_matches <- intersect(f_list_eco, basename(f.list_valRM))
    cat(
      sprintf(
        "%.2f%% (%d) of the species names of ecoregions match with names of the validation data files \n",
        100 * length(f.list_matches) / length(f.list_valRM),
        length(f.list_matches)
      )
    )
  } else {
    f.list_matches <- intersect(f_list_eco, f.list_valRM)
    cat(
      sprintf(
        "%.2f%% (%d) of the species names of ecoregions match with those in the validation data column 'sci_name'\n",
        100 * length(f.list_matches) / length(f.list_valRM),
        length(f.list_matches)
      )
    )
  }
  cat("------", "\n")
  
  df_eval <- data.frame(
    species = f.list_matches,
    nbr_pres = NA,
    nbr_pres_ecoRM = NA,
    nbr_true_pres = NA,
    nbr_false_pres = NA,
    Sen_ecoRM = NA,
    Prec_ecoRM = NA,
    type = ecoRM.dir
  )
  
  # Process each species
  for (i in seq_along(f.list_matches)) {
    if (verbose) cat(i, " Species: ", f.list_matches[i], "\n")
    
    if (is.null(mask)) {
      domain_raster <- NULL
    } else {
      domain_raster <- mask
    }
    
    # Load ecoRM raster
    ecoRM <- terra::rast(f_list_eco_full[grep(f.list_matches[i], f_list_eco_full)])
    
    # Load validation data
    if (valData.type == "TIFF") {
      valRM <- terra::rast(f.list_valRM_full[grep(f.list_matches[i], f.list_valRM_full)])
      suppressWarnings({
        if (is.na(any(all.equal(terra::crs(ecoRM), terra::crs(valRM))))) {valRM <- terra::project(valRM, ecoRM)}
      })
    } else {
      valRM <- shp_files[grep(f.list_matches[i], shp_files$sci_name), ]
      valRM <- sf::st_transform(valRM, crs = terra::crs(ecoRM))
    }
    
    ext_e <- terra::ext(ecoRM)
    ext_v <- terra::ext(valRM)
    
    if (is.null(mask)) {
      if (is.na(any(all.equal(ext_e, ext_v)))) {
        combined_extent <- terra::ext(
          min(ext_e[1], ext_v[1]),
          max(ext_e[2], ext_v[2]),
          min(ext_e[3], ext_v[3]),
          max(ext_e[4], ext_v[4])
        )
        domain_raster <- terra::extend(ecoRM, combined_extent)
        if (valData.type == "SHP") {
          valRM <- terra::rasterize(valRM, domain_raster, field = 1, background = NA)
        }
      } else {
        if (valData.type == "SHP") {
          valRM <- terra::rasterize(valRM, ecoRM, field = 1, background = NA)
        }
        domain_raster <- valRM
      }
      domain_raster[] <- 0
      domain_raster[is.na(valRM[])] <- NA
    }
    
    # Resample domain raster if resolution factor is provided
    if (!is.null(res.fact) && res.fact >= 1) {
      domain_raster <- terra::aggregate(domain_raster, fact = res.fact, fun = max, na.rm = TRUE)
    }
    
    # Ensure resolution and extent match
    ecoRM <- terra::resample(ecoRM, domain_raster, method = "near")
    valRM <- terra::resample(valRM, domain_raster, method = "near")
    
    if (!is.null(mask)) {
      ecoRM[is.na(domain_raster[])] <- NA
      valRM[is.na(domain_raster[])] <- NA
    }
    
    # Create a binarized ecoRM
    ecoRM_bin <- ecoRM
    ecoRM_bin[!is.na(domain_raster[])] <- 0
    ecoRM_bin[ecoRM == 1] <- 1
    
    # Reclassify validation raster
    if (all(range(terra::values(valRM), na.rm = TRUE) == c(0, 1))) {
      valRM_pres <- terra::classify(valRM, rbind(c(0, 0.5, 0), c(0.5, 1, 1)), include.lowest = TRUE)
    } else {
      valRM_pres <- valRM
      valRM_pres[!is.na(domain_raster[])] <- 0
      valRM_pres[valRM == 1] <- 1
    }
    
    # Burn into the overlay layer with domain_raster properties
    overlay_raster <- valRM_pres * 2 + ecoRM_bin
    overlay_raster[is.na(domain_raster[])] <- NA
    names(overlay_raster) <- f.list_matches[i]
    overlay_list[[i]] <- overlay_raster
    
    # Evaluate presence and absence
    df_eval$nbr_pres[i] <- sum(overlay_raster[] %in% c(2, 3), na.rm = TRUE)
    df_eval$nbr_pres_ecoRM[i] <- sum(overlay_raster[] %in% c(1, 3), na.rm = TRUE)
    df_eval$nbr_true_pres[i] <- sum(overlay_raster[] == 3, na.rm = TRUE)
    df_eval$nbr_false_pres[i] <- sum(overlay_raster[] == 1, na.rm = TRUE)
    
    df_eval$Sen_ecoRM[i] <- df_eval$nbr_true_pres[i] / df_eval$nbr_pres[i]
    df_eval$Prec_ecoRM[i] <- df_eval$nbr_true_pres[i] / (df_eval$nbr_true_pres[i] + df_eval$nbr_false_pres[i])
    
    # Plot the overlay raster
    if (verbose && print.map) {
      if (!dir.exists(file.path(root.dir, "Output"))) {
        dir.create(file.path(root.dir, "Output"))
      }
      
      aspect_r <- terra::nrow(overlay_raster) / terra::ncol(overlay_raster)
      colors <- c("gray", "red", "blue", "purple")
      breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
      
      pdf(file = file.path(root.dir, "Output", paste0("Evaluation_map_", f.list_matches[i], ".pdf")), height = (11 * aspect_r)+3, width = 11)
      par(mfrow = c(1, 1))
      terra::plot(
        overlay_raster,
        col = colors,
        breaks = breaks,
        legend = FALSE,
        main = paste0("Overlay validation RM and ecoRM - ", "Species: ", f.list_matches[i],"\n", 
                       "Sensitivity (TP/[TP+FA]) = ", round(df_eval$Sen_ecoRM[i], digits = 2)," & ", 
                       "Precision (TP/[TP+FP]) = ", round(df_eval$Prec_ecoRM[i], digits = 2)),
        las = 1
      )
      
      # Adding a legend
      legend("bottomright",
        legend = c(
          "Abs in both",
          "Pres in ecoRM only (FP)",
          "Pres in valRM only (FA)",
          "Pres in both (TP)"
        ),
        fill = colors,
        bg = NA,
        box.col = NA,
        inset = c(0,0.05)
      )
      dev.off()
    }
  }
  
  if (verbose) {
    cat("Cross-species mean Sen & Prec:", round(mean(rowMeans(df_eval[, c("Sen_ecoRM", "Prec_ecoRM")], na.rm = TRUE)), digits = 2), "\n")
  }
  
  if (print.map) {
    cat("### Maps have been saved to:", file.path(root.dir, "Output"), "###\n")
  }
  
  output <- list(df_eval = df_eval, overlay_list = overlay_list)
  return(output)
}
