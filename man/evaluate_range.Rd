% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_range.R
\name{evaluate_range}
\alias{evaluate_range}
\title{Evaluate the sensitivity & precision of range maps}
\usage{
evaluate_range(
  root_dir = NULL,
  valData_dir = NULL,
  ecoRM_dir = NULL,
  valData_type = NULL,
  verbose = TRUE,
  print_map = TRUE,
  mask = NULL,
  res_fact = NULL
)
}
\arguments{
\item{root_dir}{Character. Working directory to load and save target files.}

\item{valData_dir}{Numeric. Directory to validation spatial data (must have same name as data in ecoRM_dir)}

\item{ecoRM_dir}{Numeric. Directory to range maps (generated with get_range).}

\item{valData_type}{Character. Type of valData - either "SHP" or "TIFF".}

\item{verbose}{Logical (optional). Report details while running.}

\item{print_map}{Logical (optional). If verbose=TRUE should a overlap map be printed. Default is TRUE.}

\item{mask}{rast object (optional). To mask the study the focal study region. Default is TRUE.}

\item{res_fact}{Integer. Factor for coarsening the original resolution (i.e., that of the species range map)}
}
\value{
A data.frame of evaluation for all species and a list of range overlay maps. 
Precision (ppv) = number true presences (TP) / [TP + number false presences (FP)]; 
Sensitivity = number true presences (TP) / [TP + number false absences (FA)]; 
Specificity = number true absences (TA) / [TA + number false presences (FP)]; 
TSS = Sensitivity + Specificity - 1
}
\description{
Evaluates the precision (ppv), sensitivity, specificity and TSS of range maps
based on validation data, such as predictions of species distributions (SDMs)
or IUCN expert range maps. See Pinkert et al. (2023) for metrics and comparisons
of various types of range data, including expert range maps, SDMs and ecoregional
range maps. Optional functionalities include the masking of the focal study 
region (see 'mask') and aggregations of the input maps to different resolutions,
given the importance of these factors for specific applications (Pinkert et al., 2023).
}
\examples{
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
                "Sensitivity = ", round(res10km$df_eval[1,"Sen_ecoRM"], digits = 2)),
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
boxplot(Value ~ Variable + Resolution, data = long_df, 
        col = c("#FFC300", "#619CFF"), 
        names = rep(variables, 2),
        xlab = "Variable", ylab = "Value", las = 1, 
        main = "Boxplot of Sen_ecoRM and Prec_ecoRM")

# Adding legend for colors
legend("bottomright", legend = c("5km", "10km"), fill = c("#FFC300", "#619CFF"),
       title = "Resolution", bty = "n")
}
}
\references{
Pinkert, S., Sica, Y. V., Winner, K., & Jetz, W. (2023). The potential of 
ecoregional range maps for boosting taxonomic coverage in ecology and 
conservation. Ecography, 12, e06794.
}
