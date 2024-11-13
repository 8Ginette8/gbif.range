# EcoRM evaluation at different resolutions (>4min runtime)
root.dir  <- list.files(system.file(package = "gbif.range"), pattern = "extdata", full.names = TRUE)

res1km <- evaluateRM(root.dir = root.dir, 
           valData.dir = "SDM", 
           ecoRM.dir = "EcoRM",
           verbose = TRUE, 
           print.map = TRUE,
           valData.type = "TIFF", 
           mask = NULL, 
           res.fact = NULL)

res10km <- evaluateRM(root.dir = root.dir, 
                   valData.dir = "SDM", 
                   ecoRM.dir = "EcoRM",
                   verbose = TRUE, 
                   print.map = TRUE,
                   valData.type = "TIFF", 
                   mask = NULL, 
                   res.fact = 10)


# Extract and plot a specific overlay map
terra::plot(
  res10km$overlay_list[[1]],
  col = c("gray", "red", "blue", "purple"),
  breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
  legend = FALSE,
  main = paste0("Species: ", res10km$df_eval[1,1],"\n", 
                "Sensitivity = ", round(res10km$df_eval[1,6], digits = 2)," & ", 
                "Precision = ", round(res10km$df_eval[1,7], digits = 2)),
  las = 1
)

legend("bottomright",
       legend = c(
         "Abs in both",
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
res1km$df_eval$Mean_SenPrec <- (res1km$df_eval$Sen_ecoRM+res1km$df_eval$Prec_ecoRM)/2 
res10km$df_eval$Mean_SenPrec <- (res10km$df_eval$Sen_ecoRM+res10km$df_eval$Prec_ecoRM)/2 

# Combine the data frames and add a resolution column
combined_df <- bind_rows(
  res1km$df_eval %>% mutate(Resolution = "1km"),
  res10km$df_eval %>% mutate(Resolution = "10km")

)

# Gather data into long format
long_df <- combined_df %>%
  select(Sen_ecoRM, Prec_ecoRM, Mean_SenPrec, Resolution) %>%
  pivot_longer(cols = c(Sen_ecoRM, Prec_ecoRM, Mean_SenPrec), 
               names_to = "Variable", values_to = "Value")

## Plot results for different range data resolutions ##
# Sensitivity increases and precision decreases with increasing resolution
# Overall performance is slightly higher at more coarse grain
ggplot(long_df, aes(x = Variable, y = Value, fill = Resolution)) +
  geom_boxplot() +
  labs(x = "Variable", y = "Value", title = "Boxplot of Sen_ecoRM and Prec_ecoRM") +
  scale_fill_manual(values = c("10km" = "#619CFF", "1km" = "#FFC300")) +
  theme_minimal()
