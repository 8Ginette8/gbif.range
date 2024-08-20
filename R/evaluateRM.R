
library(sf)
library(raster)
library(dplyr)

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
  shp_names <-
    unique(list.files(
      file.path(root.dir, valData.dir),
      pattern = "\\.shp$",
      full.names = TRUE
    ))
  
  if (length(shp_names) != 0) {
    shp_files <- do.call(rbind, lapply(shp_names, st_read))
    if (any(duplicated(shp_files$sci_name))) {
      shp_files <- shp_files %>%
        group_by(sci_name) %>%
        summarize(geometry = st_union(geometry), .groups = "drop")
    }
    f.list_valRM <- shp_files$sci_name
  } else {
    f.list_valRM_full <-
      list.files(
        file.path(root.dir, valData.dir),
        pattern = "\\.tif$",
        full.names = TRUE,
        recursive = TRUE
      )
    f.list_valRM <- sub("\\.tif$", "", basename(f.list_valRM_full))
  }
  # Check and report match of file names
  if (length(f.list_valRM) == 0) {
    stop(
      "No names of 'valData.dir' match those of the ecoRMs.
         Assuming that your species names of the ecoRMs are separated by a space and followed by '.tif'.
         Text before the species names (all separated with underscores) are ignored. \n"
    )
  }
  
  cat("-Note-","\n")
  if (valData.type == "TIF") {
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
  cat("---","\n")
  
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
    if (verbose)
      cat(i, " Species: ", f.list_matches[i], "\n")
    
    if (is.null(mask)) {
      domain_raster <- NULL
    } else {
      domain_raster <- mask
    }
    
    # Load ecoRM raster
    ecoRM <-
      raster(f_list_eco_full[grep(f.list_matches[i], f_list_eco_full)])
    
    # Load validation data
    if (valData.type == "TIF") {
      valRM <-
        raster(f.list_valRM_full[grep(f.list_matches[i], f.list_valRM_full)])
      if (is.na(any(all.equal(crs(ecoRM), crs(valRM)))))
        valRM <- projectRaster(valRM, ecoRM)
    } else {
      valRM <- shp_files[grep(f.list_matches[i], shp_files$sci_name),]
      valRM <- as(st_transform(valRM, crs(ecoRM)), "Spatial")
    }
    
    
    if (is.null(mask)) {
      if (is.na(any(all.equal(extent(ecoRM), extent(valRM))))) {
        combined_extent <- extent(
          min(extent(ecoRM)@xmin, extent(valRM)@xmin),
          max(extent(ecoRM)@xmax, extent(valRM)@xmax),
          min(extent(ecoRM)@ymin, extent(valRM)@ymin),
          max(extent(ecoRM)@ymax, extent(valRM)@ymax)
        )
        domain_raster <- extend(ecoRM, combined_extent)
        if (valData.type == "SHP") {
          valRM <- raster::rasterize(valRM,
                                     domain_raster,
                                     field = 1,
                                     background = NA)
        }
      } else  {
        if (valData.type == "SHP") {
          valRM <- raster::rasterize(valRM,
                                     ecoRM,
                                     field = 1,
                                     background = NA)
        }
        domain_raster  <- valRM
      }
      domain_raster[] <- 0
      domain_raster[is.na(valRM[])] <- NA
    }
    
    # Resample domain raster if resolution factor is provided
    if (!is.null(res.fact) | res.fact >= 1) {
      domain_raster <-
        aggregate(
          domain_raster,
          fact = res.fact,
          FUN = function(x)
            max(x, na.rm = TRUE)
        )
    }
    
    # Ensure resolution and extent match
    ecoRM <- raster::resample(ecoRM, domain_raster, method = "ngb")
    valRM <- raster::resample(valRM, domain_raster, method = "ngb")
    
    if (!is.null(mask)) {
      ecoRM[is.na(domain_raster[])] <- NA
      valRM[is.na(domain_raster[])] <- NA
    }
    
    # Create a binarized ecoRM
    ecoRM_bin <- ecoRM
    ecoRM_bin[!is.na(domain_raster[])] <- 0
    ecoRM_bin[ecoRM == 1] <- 1
    
    # Reclassify validation raster
    if (all(range(getValues(valRM), na.rm = TRUE) == c(0, 1))) {
      valRM_pres <-  raster::reclassify(valRM, c(0, 0.5, 0, 0.5, 1, 1))
    }
    
    if (!all(range(getValues(valRM), na.rm = TRUE) == c(0, 1))) {
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
    df_eval$nbr_pres[i] <-
      sum(overlay_raster[] %in% c(2, 3), na.rm = TRUE)
    df_eval$nbr_pres_ecoRM[i] <-
      sum(overlay_raster[] %in% c(1, 3), na.rm = TRUE)
    df_eval$nbr_true_pres[i] <-
      sum(overlay_raster[] == 3, na.rm = TRUE)
    df_eval$nbr_false_pres[i] <-
      sum(overlay_raster[] == 1, na.rm = TRUE)
    
    df_eval$Sen_ecoRM[i] <-
      df_eval$nbr_true_pres[i] / df_eval$nbr_pres[i]
    df_eval$Prec_ecoRM[i] <-
      df_eval$nbr_true_pres[i] / (df_eval$nbr_true_pres[i] + df_eval$nbr_false_pres[i])
    
    # Plot the overlay raster
    if (verbose && print.map) {
      if (!dir.exists(file.path(root.dir, "Output"))) {
        dir.create(file.path(root.dir, "Output"))
      }
      
      aspect_r <- overlay_raster@nrows / overlay_raster@ncols
      colors <- c("gray", "red", "blue", "purple")
      breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
      
      pdf(file=file.path(root.dir, "Output" ,paste0("Evaluation_map_", f.list_matches[i],".pdf")), height=(11*aspect_r)+4, width=11)
      par(mfrow=c(1,1))
      plot(
        overlay_raster,
        col = colors,
        breaks = breaks,
        legend = FALSE,
        main = sprintf("Overlay validation RM and ecoRM\n%s", f.list_matches[i]),
        las = 1
      )
      
      # Adding a legend
      legend(
        "bottomright",
        legend = c(
          "Abs in both",
          "Pres in ecoRM only (FP)",
          "Pres in valRM only (FA)",
          "Pres in both (TP)"
        ),
        fill = colors,
        bg = NA,
        box.col = NA
      )
      legend(
        "topleft",
        legend = c(
          sprintf("Sensitivity (TP/[TP+FA]) = %.2f", df_eval$Sen_ecoRM[i]),
          sprintf("Precision (TP/[TP+FP]) = %.2f", df_eval$Prec_ecoRM[i])
        ),
        bg = NA,
        box.col = NA
      )
      dev.off()
    }
  }
  
  if (verbose) {
    cat("Cross-species mean Sen & Prec:",
        round(mean(rowMeans(
          df_eval[, c("Sen_ecoRM", "Prec_ecoRM")], na.rm = TRUE
        )), digits = 2),
        "\n")
  }
  
  if (print.map) {
    cat("### Maps have been saved to:",
        file.path(root.dir, "Output"),
        "###\n")
  }
  output <- list(df_eval = df_eval, overlay_list = overlay_list)
  return(output)
}