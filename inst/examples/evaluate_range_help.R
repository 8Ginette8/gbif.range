\dontrun{
# EcoRM evaluation at different resolutions (>4min runtime)
root.dir  <- list.files(system.file(package = "gbif.range"), pattern = "extdata", full.names = TRUE)

res5km <- evaluate_range(root_dir = root.dir, 
           valData_dir = "SDM", 
           ecoRM_dir = "EcoRM",
           verbose = TRUE, 
           print_map = TRUE,
           valData_type = "TIFF", 
           mask = NULL, 
           res_fact = NULL)

res10km <- evaluate_range(root_dir = root.dir, 
                   valData_dir = "SDM", 
                   ecoRM_dir = "EcoRM",
                   verbose = TRUE, 
                   print_map = TRUE,
                   valData_type = "TIFF", 
                   mask = NULL, 
                   res_fact = 2)


# Extract and plot a specific overlay map
terra::plot(
  res10km$overlay_list[[1]],
  col = c("gray", "red", "blue", "purple"),
  breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
  legend = FALSE,
  main = paste0("Species: ", res10km$df_eval[1,1],"\n", 
                "Precision = ", round(res10km$df_eval[1,"Prec_ecoRM"], digits = 2)," & ", 
                "TSS = ", round(res10km$df_eval[1,"TSS_ecoRM"], digits = 2)),
  las = 1
)

legend("bottomright",
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
res5km$df_eval$Mean_SenPrec <- (res5km$df_eval$Sen_ecoRM + res5km$df_eval$Prec_ecoRM) / 2
res10km$df_eval$Mean_SenPrec <- (res10km$df_eval$Sen_ecoRM + res10km$df_eval$Prec_ecoRM) / 2

# Combine the data frames and add a Resolution column
combined_df <- rbind(
  cbind(res5km$df_eval, Resolution = "1km"),
  cbind(res10km$df_eval, Resolution = "10km")
)

# Convert to long format
variables <- c("Sen_ecoRM", "Prec_ecoRM", "Mean_SenPrec")
long_df <- data.frame(
  Variable = rep(variables, each = nrow(combined_df)),
  Value = unlist(combined_df[variables]),
  Resolution = rep(combined_df$Resolution, times = length(variables))
)

# Plot boxplots using base R
boxplot(Value ~ Variable + Resolution, data = long_df, 
        col = c("#FFC300", "#619CFF"), 
        names = rep(variables, 2),
        xlab = "Variable", ylab = "Value", las = 1, 
        main = "Boxplot of Sen_ecoRM and Prec_ecoRM")

# Adding legend for colors
legend("bottomright", legend = c("1km", "10km"), fill = c("#FFC300", "#619CFF"),
       title = "Resolution", bty = "n")
}