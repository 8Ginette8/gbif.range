\dontrun{
###########################################
### Example plot
###########################################

## EcoRM evaluation at different resolutions (>4min runtime)
root.dir  <- list.files(
    system.file(package = "gbif.range"),
    pattern = "extdata",
    full.names = TRUE
)

res5km <- evaluate_range(
    root_dir = root.dir, 
    valData_dir = "SDM", 
    ecoRM_dir = "EcoRM",
    verbose = TRUE, 
    print_map = TRUE,
    valData_type = "TIFF", 
    mask = NULL, 
    res_fact = NULL
)

res10km <- evaluate_range(
    root_dir = root.dir, 
    valData_dir = "SDM", 
    ecoRM_dir = "EcoRM",
    verbose = TRUE, 
    print_map = TRUE,
    valData_type = "TIFF", 
    mask = NULL, 
    res_fact = 2
)


## Extract and plot a specific overlay map
terra::plot(
    res10km$overlay_list[[1]],
    col = c("gray", "red", "blue", "purple"),
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
    legend = FALSE,
    main = paste0("Species: ",
        res10km$df_eval[1, 1],"\n", 
        "Precision = ",
        round(res10km$df_eval[1, "Prec_ecoRM"], digits = 2)," & ", 
          "Sensitivity = ",
        round(res10km$df_eval[1, "Sen_ecoRM"], digits = 2)),
    las = 1
)

legend(
    "bottomright",
    legend = c(
      "Abs in both (TA)",
      "Pres in ecoRM only (FP)",
      "Pres in valRM only (FA)",
      "Pres in both (TP)"
    ),
    fill = c("gray", "red", "blue", "purple"),
    bg = NA,
    box.col = NA,
    inset = c(0,0.2)
)

## Compare sensitivity and precision outputs ##
  # Calculate overall performance (e.g. mean sensitivity & precision)
res5km$df_eval$Mean_SenPrec <-
    (res5km$df_eval$Sen_ecoRM + res5km$df_eval$Prec_ecoRM) / 2
res10km$df_eval$Mean_SenPrec <-
    (res10km$df_eval$Sen_ecoRM + res10km$df_eval$Prec_ecoRM) / 2

  # Combine the data frames and add a Resolution column
combined.df <- rbind(
    cbind(res5km$df_eval, Resolution = "5km"),
    cbind(res10km$df_eval, Resolution = "10km")
)

  # Convert to long format
variables <- c("Sen_ecoRM", "Prec_ecoRM", "Mean_SenPrec")
long_df <- data.frame(
    Variable = rep(variables, each = nrow(combined.df)),
    Value = unlist(combined.df[variables]),
    Resolution = rep(combined.df$Resolution, times = length(variables))
)

  # Plot boxplots using base R
boxplot(
    Value ~ Variable + Resolution,
    data = long_df, 
    col = c("#FFC300", "#619CFF"), 
    names = rep(variables, 2),
    xlab = "Variable",
    ylab = "Value",
    las = 1, 
    main = "Boxplot of Sen_ecoRM and Prec_ecoRM"
)

  # Adding legend for colors
legend(
    "bottomright",
    legend = c("5km", "10km"),
    fill = c("#FFC300", "#619CFF"),
    title = "Resolution",
    bty = "n"
)

###########################################
### Package manuscript plot (Fig 2c-d)
###########################################

# Root and package
root_dir <- list.files(
    system.file(package = "gbif.range"),
    pattern = "extdata",
    full.names = TRUE
)
if (!dir.exists(file.path(root_dir, "fig_plots"))) {
    dir.create(file.path(root_dir, "fig_plots"))
}
if (!requireNamespace("colorspace", quietly = TRUE)) {
  install.packages("colorspace")
}

###########
##### Plant #########
###########

# Preliminary
shp.lonlat <- terra::vect(
  paste0(system.file(package = "gbif.range"),"/extdata/shp_lonlat.shp")
)
obs.arcto <- get_gbif(
    "Arctostaphylos alpinus",
    geo = shp.lonlat,
    grain = 1
)
rst <- terra::rast(
  paste0(system.file(package = "gbif.range"),"/extdata/rst.tif")
)
my.eco <- make_ecoregion(rst, 200)
range.arcto <- get_range(
    occ_coord = obs.arcto,
    bioreg = my.eco,
    bioreg_name = "EcoRegion",
    res = 0.05
)
ext.temp <- terra::ext(range.arcto$rangeOutput)
ext.temp <- c(ext.temp[1]-0.6, ext.temp[2]+0.05,
                ext.temp[3]-0.05, ext.temp[4]+0.05)
spdf.world <- terra::vect(
  rnaturalearth::ne_countries(scale = 10,returnclass = "sf")
)
world.local <- terra::crop(spdf.world,ext.temp)
world.local.ar <- terra::aggregate(world.local)
r.arcto <- terra::mask(range.arcto$rangeOutput, world.local.ar)

# Plot plant
png(
   paste0(
        root_dir,
        "/fig_plots/fig2_arcto_ind.png"
  ),
  width = 100,
  height = 70,
  unit = "cm",
  res = 100,
  pointsize = 110
)
par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 1, cex = 0.5)
terra::plot(world.local.ar, col = "#dce8dc", axes = FALSE, lwd = 2)
terra::plot(
  terra::as.polygons(
    r.arcto
  ),
  lwd = 0.1,
  col = "#63636350",
  add = TRUE
)
terra::plot(
  res5km$overlay_list[[1]],
  col = c("#d1c845", "#e39c59", "#8d60cc", "#6cba6c"),
  breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
  axes=FALSE,
  legend = FALSE,
  las = 1,
  add = TRUE
)
text(
  6.6,43.2,
  paste("Mean TSS =", round(res5km$df_eval[1,"TSS_ecoRM"],2)),
  cex = 1.5,
  font = 2
)
text(
  7.2,42.8,
  paste("Mean Precision =", round(res5km$df_eval[1,"Prec_ecoRM"],2)),
  cex = 1.5,
  font = 2
)
terra::plot(
  world.local.ar,
  border = "#383d38",
  axes = FALSE,
  lwd = 2,
  add = TRUE
)
dev.off()

###########
##### Tiger #########
###########

# Preliminary
obs.pt <- get_gbif(sp_name = "Panthera tigris")
eco.terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)
range.tiger <- get_range(
    occ_coord = obs.pt,
    bioreg = eco.terra,
    bioreg_name = "ECO_NAME"
)
ext.temp <- terra::ext(range.tiger$rangeOutput)
ext.temp <- c(ext.temp[1]-0.2, ext.temp[2]+0.2,
                ext.temp[3]-0.2, ext.temp[4]+0.2)
world.local <- terra::crop(spdf.world, ext.temp)
world.local.ti <- terra::aggregate(world.local)

# Plot tiger
png(
  paste0(
        root_dir,
        "/fig_plots/fig2_tiger_ind.png"
  ),
  width = 100,
  height = 70,
  unit = "cm",
  res = 100,
  pointsize = 110
)
par(mfrow = c(1,1), mar = c(5,5,5,20), lwd = 1, cex = 0.5)
terra::plot(world.local.ti, col = "#dce8dc", axes = FALSE, lwd = 2)
toPlot = terra::mask(res5km$overlay_list[[6]], world.local.ti)
terra::plot(
  toPlot,
  col = c("#d1c845", "#e39c59", "#8d60cc", "#6cba6c"),
  breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
  axes = FALSE,
  legend = FALSE,
  las = 1,
  add = TRUE
)
text(
  86,49,
  paste("Mean TSS =",round(res5km$df_eval[6,"TSS_ecoRM"],2)),
  cex = 1.5,
  font = 2
)
text(
  90.4,45.4,
  paste("Mean Precision =",round(res5km$df_eval[6,"Prec_ecoRM"],2)),
  cex = 1.5,
  font = 2
)
terra::plot(world.local.ti,
  border = "#383d38",
  axes = FALSE,
  lwd = 2,
  add = TRUE
)
legend(
  "bottomleft",
  legend = c(
    "True Presences",
    "True Absences",
    "False Presences",
    "False Absences"
  ),
  fill = c("#6cba6c", "#d1c845", "#e39c59", "#8d60cc"),
  bg = NA,
  box.col = NA,
  inset = c(0.13,0),
  cex = 1.1,
  x.intersp = 0.2
)
dev.off()

}