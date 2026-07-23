# Evaluate Range Maps Against Independent Validation Data

Compare range maps produced by
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
with external validation data, such as species distribution models
(SDMs) or expert-derived range maps, and summarize precision,
sensitivity, specificity, and TSS.

## Usage

``` r
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
```

## Arguments

- root_dir:

  Character. String giving the root directory that contains both the
  generated range maps and the validation data.

- valData_dir:

  Character. String giving the directory, relative to `root_dir`,
  containing the validation data.

- ecoRM_dir:

  Character. String giving the directory, relative to `root_dir`,
  containing range maps generated with
  [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md).

- valData_type:

  Character. String indicating the expected validation-data format:
  `"SHP"` or `"TIFF"`.

- verbose:

  Logical. Should progress information be printed while the comparison
  is running?

- print_map:

  Logical. If `TRUE`, write a PDF overlay map for each evaluated
  species.

- mask:

  Optional `SpatRaster`. Used as a study-area mask and common comparison
  domain.

- res_fact:

  Optional integer. Aggregation factor used to coarsen the native
  resolution before comparison.

## Value

A list with two elements: `df_eval`, a data frame containing per-species
evaluation statistics, and `overlay_list`, a list of raster overlays
used for plotting and inspection.

## Details

TIFF validation files must have file names matching the species names of
the range maps. Shapefile-based validation data must include a column
named `sci_name` with matching species names.

The function can optionally mask the focal study region and aggregate
maps to coarser resolutions before calculating evaluation metrics, which
is useful when comparing products with different native resolutions.

## References

Pinkert, S., Sica, Y. V., Winner, K., & Jetz, W. (2023). The potential
of ecoregional range maps for boosting taxonomic coverage in ecology and
conservation. Ecography, 12, e06794.
[doi:10.1111/ecog.06794](https://doi.org/10.1111/ecog.06794)

## See also

[`get_range`](https://8ginette8.github.io/gbif.range/reference/get_range.md)()
to generate the range maps being evaluated, and
[`cv_range`](https://8ginette8.github.io/gbif.range/reference/cv_range.md)()
for cross-validation against the original occurrence data instead of
external validation data.

## Examples

``` r
## EcoRM evaluation at different resolutions (<10sec runtime)
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
    print_map = FALSE,
    valData_type = "TIFF", 
    mask = NULL, 
    res_fact = NULL
)
#> -Note-
#> 100.00% (12) of the species names of ecoregions match with names of the validation data files
#> ------
#> 1 Species: Arctostaphylos_alpinus
#> 2 Species: Atropa_bella-donna
#> 3 Species: Centaurea_montana
#> 4 Species: Cirsium_rivulare
#> 5 Species: Cupressus_sempervirens
#> 6 Species: Panthera_tigris
#> 7 Species: Pinus_nigra
#> 8 Species: Plantago_alpina
#> 9 Species: Plantago_media
#> 10 Species: Prunella_vulgaris
#> 11 Species: Salvia_pratensis
#> 12 Species: Tragopogon_dubius
#> Cross-species mean Prec & Sensitivity: 0.75

res10km <- evaluate_range(
    root_dir = root.dir, 
    valData_dir = "SDM", 
    ecoRM_dir = "EcoRM",
    verbose = TRUE, 
    print_map = FALSE,
    valData_type = "TIFF", 
    mask = NULL, 
    res_fact = 2
)
#> -Note-
#> 100.00% (12) of the species names of ecoregions match with names of the validation data files
#> ------
#> 1 Species: Arctostaphylos_alpinus
#> 2 Species: Atropa_bella-donna
#> 3 Species: Centaurea_montana
#> 4 Species: Cirsium_rivulare
#> 5 Species: Cupressus_sempervirens
#> 6 Species: Panthera_tigris
#> 7 Species: Pinus_nigra
#> 8 Species: Plantago_alpina
#> 9 Species: Plantago_media
#> 10 Species: Prunella_vulgaris
#> 11 Species: Salvia_pratensis
#> 12 Species: Tragopogon_dubius
#> Cross-species mean Prec & Sensitivity: 0.75


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

graphics::legend(
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
graphics::boxplot(
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
graphics::legend(
    "bottomright",
    legend = c("5km", "10km"),
    fill = c("#FFC300", "#619CFF"),
    title = "Resolution",
    bty = "n"
)
```
