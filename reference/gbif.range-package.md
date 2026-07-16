# gbif.range: Tools for GBIF Retrieval and Ecoregion-Based Species Range Mapping

`gbif.range` provides a workflow to retrieve occurrence records from
GBIF, clean and filter them for spatial analyses, assemble or download
ecoregion layers, generate ecologically informed species range maps, and
evaluate the resulting maps against independent data.

## Details

The package is designed around a typical workflow: (1) inspect or
download GBIF records with
[`get_gbif_count()`](https://8ginette8.github.io/gbif.range/reference/get_gbif_count.md)
and
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md),
(2) retrieve taxonomic information with
[`get_status()`](https://8ginette8.github.io/gbif.range/reference/get_status.md),
(3) load packaged ecoregions with
[`read_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)
or create custom ones with
[`make_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md),
(4) build range maps with
[`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md),
and (5) evaluate those maps with
[`evaluate_range()`](https://8ginette8.github.io/gbif.range/reference/evaluate_range.md)
or
[`cv_range()`](https://8ginette8.github.io/gbif.range/reference/cv_range.md).

Additional helpers include
[`get_doi()`](https://8ginette8.github.io/gbif.range/reference/get_doi.md)
for creating GBIF-derived dataset DOIs,
[`obs_filter()`](https://8ginette8.github.io/gbif.range/reference/obs_filter.md)
for grid-based thinning,
[`make_tiles()`](https://8ginette8.github.io/gbif.range/reference/make_tiles.md)
for splitting study extents into GBIF-ready polygons, and a disk-based
workflow built around
[`split_gbif_by_species()`](https://8ginette8.github.io/gbif.range/reference/split_gbif_by_species.md),
[`species_csvs_to_ranges()`](https://8ginette8.github.io/gbif.range/reference/species_csvs_to_ranges.md),
and
[`read_range_rds()`](https://8ginette8.github.io/gbif.range/reference/read_range_rds.md)
for very large downloaded GBIF tables.

## See also

Useful links:

- <https://github.com/8Ginette8/gbif.range>

- <https://8ginette8.github.io/gbif.range/>

- Report bugs at <https://github.com/8Ginette8/gbif.range/issues>

## Author

**Maintainer**: Yohann Chauvier <yohann.chauvier@wsl.ch>
([ORCID](https://orcid.org/0000-0001-9399-3192))

Authors:

- Yohann Chauvier <yohann.chauvier@wsl.ch>
  ([ORCID](https://orcid.org/0000-0001-9399-3192))

- Oskar Hagen <oskar@hagen.bio>
  ([ORCID](https://orcid.org/0000-0002-7931-6571))

- Stefan Pinkert ([ORCID](https://orcid.org/0000-0002-8348-2337))

- Camille Albouy ([ORCID](https://orcid.org/0000-0003-1629-2389))

- Fabian Fopp ([ORCID](https://orcid.org/0000-0003-0648-8484))

- Philipp Brun ([ORCID](https://orcid.org/0000-0002-2750-9793))

- Patrice Descombes ([ORCID](https://orcid.org/0000-0002-3760-9907))

- Florian Altermatt ([ORCID](https://orcid.org/0000-0002-4831-6958))

- Loïc Pellissier ([ORCID](https://orcid.org/0000-0002-2289-8259))

- Katalin Csillery ([ORCID](https://orcid.org/0000-0003-0039-9296))

## Examples

``` r
# -------------------------------------------------------------------------
# 1. Minimal in-memory workflow with custom ecoregions
# -------------------------------------------------------------------------

# Create two simple environmental layers on a small synthetic study area.
r1 <- terra::rast(ncols = 20, nrows = 20, xmin = 0, xmax = 10,
  ymin = 0, ymax = 10)
terra::values(r1) <- rep(seq(0, 1, length.out = 20), each = 20)
r2 <- terra::rast(r1)
terra::values(r2) <- rep(seq(0, 1, length.out = 20), times = 20)
env <- c(r1, r2)

# Derive custom ecoregions and define a small occurrence table in memory.
eco <- make_ecoreg(env = env, nclass = 4)
#> CLARA algorithm processing... 
#> Generating polygons... 
occ <- data.frame(
  decimalLongitude = c(0.5, 1.2, 2.4, 3.6, 2.8, 6.0, 7.2, 7.4, 7.6, 8.8),
  decimalLatitude = c(1.0, 0.1, 2.3, 2.5, 2.7, 5.0, 7.1, 7.3, 6.5, 7.7)
)

# Build the range directly from the occurrence table and ecoregions.
range_obj <- get_range(
  occ_coord = occ,
  ecoreg = eco,
  verbose = FALSE
)

# Plot the predicted range and overlay the occurrence points.
terra::plot(range_obj$rangeOutput, col = 3, main = "Range Map")
graphics::points(occ, pch = 4)


# \donttest{
# -------------------------------------------------------------------------
# 2. Typical online workflow with GBIF data
# -------------------------------------------------------------------------

# Download GBIF occurrences for one species.
obs <- get_gbif("Panthera tigris", grain = 25)
#> |--------------------------------------------|
#> | Total number (all records)    :       8007 |
#> | Kept records                  :       5457 |
#> |--------------------------------------------|
#> | Kept records according to parameters:
#> | spatial_issue = FALSE, has_xy = TRUE
#> 
#> ...GBIF records of Panthera tigris: download starting...
#> ------------- #1 (100%..)               
#> 
#> ...Records (XY) filtering summary:
#> -----------------------------------------------
#>                     step removed remaining
#>          Grain filtering    4082      1375
#>       Duplicated records     793       582
#>          Absence records       0       582
#>          Basis selection      72       510
#>  Establishment selection       0       510
#>               Time frame       0       510
#>        Identical records       0       510
#>         Raster centroids       0       510
#> 
#> Initial records         : 5457
#> Total removed           : 4947
#> Final records (XY)      : 510
#> -----------------------------------------------
#> Final records (no XY)   : 0

# Inspect the GBIF backbone interpretation used by the package.
status <- get_status("Panthera tigris")
status
#>           canonicalName    rank gbif_key                   scientificName
#> 5219416 Panthera tigris SPECIES  5219416 Panthera tigris (Linnaeus, 1758)
#> 4969819    Felis tigris SPECIES  4969819      Felis tigris Linnaeus, 1758
#>         gbif_status    Genus  Family     Order    Class   Phylum IUCN_status
#> 5219416    ACCEPTED Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#> 4969819     SYNONYM Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#>         sp_nameMatch
#> 5219416        INPUT
#> 4969819        EXACT

# Load a packaged terrestrial ecoregion layer and build the range.
eco_terra <- read_ecoreg("eco_terra")
tiger_range <- get_range(
  occ_coord = obs,
  ecoreg = eco_terra,
  ecoreg_name = "ECO_NAME"
)
#> ## Start of computation for species:  Panthera tigris  ### 
#> 27 outlier's from 508 | proportion from total points: 5%
#> ecoregion 1  of  27 :  Brahmaputra Valley Semi-Evergreen Forests 
#> ecoregion 2  of  27 :  Central Deccan Plateau Dry Deciduous Forests 
#> ecoregion 3  of  27 :  Chao Phraya Freshwater Swamp Forests 
#> ecoregion 4  of  27 :  Eastern Highlands Moist Deciduous Forests 
#> ecoregion 5  of  27 :  Eastern Himalayan Alpine Shrub And Meadows 
#> ecoregion 6  of  27 :  Eastern Himalayan Broadleaf Forests 
#> ecoregion 7  of  27 :  Eastern Java-Bali Rain Forests 
#> ecoregion 8  of  27 :  Himalayan Subtropical Broadleaf Forests 
#> ecoregion 9  of  27 :  Himalayan Subtropical Pine Forests 
#> ecoregion 10  of  27 :  Kayah-Karen Montane Rain Forests 
#> ecoregion 11  of  27 :  Khathiar-Gir Dry Deciduous Forests 
#> ecoregion 12  of  27 :  Luang Prabang Montane Rain Forests 
#> ecoregion 13  of  27 :  Meghalaya Subtropical Forests 
#> ecoregion 14  of  27 :  Mizoram-Manipur-Kachin Rain Forests 
#> ecoregion 15  of  27 :  Narmada Valley Dry Deciduous Forests 
#> ecoregion 16  of  27 :  Orissa Semi-Evergreen Forests 
#> ecoregion 17  of  27 :  Peninsular Malaysian Montane Rain Forests 
#> ecoregion 18  of  27 :  Peninsular Malaysian Rain Forests 
#> ecoregion 19  of  27 :  South Deccan Plateau Dry Deciduous Forests 
#> ecoregion 20  of  27 :  South Western Ghats Moist Deciduous Forests 
#> ecoregion 21  of  27 :  South Western Ghats Montane Rain Forests 
#> ecoregion 22  of  27 :  Sumatran Lowland Rain Forests 
#> ecoregion 23  of  27 :  Sumatran Montane Rain Forests 
#> ecoregion 24  of  27 :  Sundarbans Mangroves 
#> ecoregion 25  of  27 :  Tenasserim-South Thailand Semi-Evergreen Rain Forests 
#> ecoregion 26  of  27 :  Terai-Duar Savanna And Grasslands 
#> ecoregion 27  of  27 :  Upper Gangetic Plains Moist Deciduous Forests 
#> ## End of computation for species:  Panthera tigris  ### 

# Plot the predicted terrestrial range and the GBIF occurrences.
terra::plot(
  tiger_range$rangeOutput,
  col = 3,
  main = paste("Range:", obs$scientificName[1]))

graphics::points(
  obs$decimalLongitude,
  obs$decimalLatitude,
  pch = 4,
  col = grDevices::rgb(1, 0, 1, 0.2)
)


# }

# \donttest{
# -------------------------------------------------------------------------
# 3. Large downloaded GBIF table already stored on disk
# -------------------------------------------------------------------------

if (requireNamespace("data.table", quietly = TRUE)) {
  # Use the bundled GBIF-style example file as a stand-in for a large download.
  gbif_file <- system.file("extdata", "occ_example_2sps.csv", package = "gbif.range")

  # Keep each stage in its own temporary folder so outputs are easy to inspect.
  split_dir <- file.path(tempdir(), "gbif_pkg_split")
  occ_dir <- file.path(tempdir(), "gbif_pkg_occ_min")
  range_dir <- file.path(tempdir(), "gbif_pkg_ranges")

  # Remove earlier temporary outputs so rerunning the example starts cleanly.
  unlink(split_dir, recursive = TRUE)
  unlink(occ_dir, recursive = TRUE)
  unlink(range_dir, recursive = TRUE)

  # Split the input table into one occurrence file per speciesKey without
  # loading the full file into memory.
  split_summary <- split_gbif_by_species(
    input_file = gbif_file,
    outdir = split_dir,
    chunk_size = 10,
    sep_in = "\t",
    sep_out = "\t",
    overwrite = TRUE,
    verbose = FALSE
  )

  split_summary[, c("species_name", "n_records", "species_file")]

  # Use the packaged terrestrial ecoregions for the batch range step.
  # species_csvs_to_ranges() will resolve "eco_terra" with read_ecoreg().
  range_summary <- species_csvs_to_ranges(
    species_dir = split_dir,
    ecoreg = "eco_terra",
    ecoreg_name = "ECO_NAME",
    outdir = range_dir,
    occ_outdir = occ_dir,
    occ_save_as = "tsv",
    range_save_as = "rds",
    sep_in = "\t",
    overwrite = TRUE,
    degrees_outlier = 30,
    clust_pts_outlier = 2,
    buff_width_point = 1,
    buff_incrmt_pts_line = 0.1,
    buff_width_polygon = 1,
    format = "SpatVector",
    verbose = FALSE
  )

  range_summary[, c("species_name", "n_points", "range_file")]

  # Read the saved ranges back from disk and plot both species together.
  range_one <- read_range_rds(range_summary$range_file[1])
  range_two <- read_range_rds(range_summary$range_file[2])

  combined_ext <- terra::ext(
    min(terra::xmin(range_one$rangeOutput), terra::xmin(range_two$rangeOutput)),
    max(terra::xmax(range_one$rangeOutput), terra::xmax(range_two$rangeOutput)),
    min(terra::ymin(range_one$rangeOutput), terra::ymin(range_two$rangeOutput)),
    max(terra::ymax(range_one$rangeOutput), terra::ymax(range_two$rangeOutput))
  )

  terra::plot(combined_ext, col = NA, legend = FALSE)
  terra::plot(range_one$rangeOutput, col = grDevices::rgb(0.1, 0.6, 0.2, 0.5), add = TRUE)
  terra::plot(range_two$rangeOutput, col = grDevices::rgb(0.8, 0.3, 0.1, 0.5), add = TRUE)
}


# }
```
