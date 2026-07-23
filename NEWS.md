# gbif.range 1.8.0
* Move manuscript plots in a dedicated vignette
* Updated the examples and vignettes
* Lower the size and time length of the examples
* Debugging
* Documentation updated

# gbif.range 1.7.1
* Updated documentation

# gbif.range 1.7.0
* Updated documentation
* Updated package framework
* Package structure and verbose updated

# gbif.range 1.6.3

* Added pkgdown website.
* Added Part 0 vignette "Getting Started" with pre-computed figures.
* Updated all vignettes with pre-computed figures.
* Updated `README.md`.
* Updated parameters in `get_gbif()` and `get_status()`.

# gbif.range 1.6.2

* Updated examples and documentation.
* Updated vignettes.
* Updated `get_status()`.
* Corrected example issues and code incompatibility.

# gbif.range 1.6.1

* Updated `get_status()`: new `level` parameter.
* Minor naming fix.
* Updated `get_gbif()`: new progress bars.

# gbif.range 1.6.0

* Added disk-based batch workflow for large multi-species GBIF exports:
  `split_gbif_by_species()`, `species_csvs_to_ranges()`, `read_range_rds()`.
* Added `Collate` field to `DESCRIPTION` for explicit R file load ordering.
* Added three focused workflow vignettes.
* Updated `README.md` with Vignettes section.

# gbif.range 1.5.3

* Improved documentation of package and `get_status()`.
* Updated `get_gbif()` with new backend support (occ_download parameters)

# gbif.range 1.5.2

* Added `area_data`: bundled dataset of `gbif.range`- vs. IUCN-derived range
  area estimates, added alongside new examples mirroring the draft paper plots.
* Updated test routine.
* Updated `get_status()`.

# gbif.range 1.5.1

* Clarified taxonomic harmonization in `get_gbif()` documentation.
* Updated `gbif.range` Rd.

# gbif.range 1.5.0

* Polished documentation and CI.
* Updated evaluation functions and examples.
* Added `get_gbif_count()`: estimate record volume before downloading.
* Changing functions name to: make_ecoreg(), get_ecoreg(),
  check_and_get_ecoreg()

# gbif.range 1.4.7

* Added package-level help page.
* Updated `get_gbif()`.
* Used `&&` in GBIF status checks.
* Updated `.gitignore` and `README.md`.

# gbif.range 1.4.0

* Added `evaluate_range()`: validate range maps against independent distribution
  data (SDMs, IUCN polygons) with precision, sensitivity, specificity, and TSS.
* Added `cv_range()`: cross-validate a `get_range()` output against its own
  occurrence data using spatial or random folds.
* Added `make_blocks()`: split observations into balanced random or spatially
  structured folds for cross-validation workflows.
* Added `area_data`: bundled dataset of `gbif.range`- vs. IUCN-derived range
  area estimates for validation examples.
* Added internal helper functions (`helpers.R`) for argument checking and
  shared utilities across functions.
* Introduced formal R5 reference classes `getRange` and `getGBIF` (`classes.R`)
  to store function outputs with their original arguments.
* Added `check_and_get_bioreg()` and `get_bioreg()` as helpers for ecoregion
  download and caching.

# gbif.range 1.1.0

* Added `make_ecoregion()`: build custom ecoregion layers from environmental
  rasters via k-means clustering.
* Moved dependencies from `Depends` to `Imports` for cleaner namespace handling.
* Added `sf` and `cluster` as dependencies.
* Expanded `get_range()` with additional ecoregion flexibility and resolution
  control via the `res` argument.

# gbif.range 1.0.0

* First stable release under the name `gbif.range` on GitHub.
* Full migration from `raster` to `terra` (SpatRaster/SpatVector compatibility).
* Renamed `get_taxonomy()` to `get_status()`: added IUCN Red List status
  retrieval and infra-specific taxa (subspecies, varieties) lookup.
* Improved `get_gbif()` synonym handling and tiling robustness.
* Added `read_ecoreg()` and `ecoreg_list` for bundled ecoregion management.

# gbif.range 0.2.0

* Renamed package from `wsl.gbif` to `gbif.range`.
* Same for function containing 'wsl': wsl_doi() to get_doi(), wsl_gbif() to
  get_gbif, wsl_obs_filter() to obs_filter()
* Added `get_range()`: ecoregion-constrained species range inference.
* Added `conv_function()`: internal polygon builder used by `get_range()`.
* Added `get_taxonomy()` (later renamed `get_status()`): GBIF backbone
  taxonomy inspection including accepted names and synonyms.
* Expanded `get_gbif()` with dynamic moving-window tiling for > 100,000
  records and improved synonym-aware downloads.
* Added `ClusterR`, `FNN`, `geometry`, `mclust`, and `rgeos` as dependencies.

# wsl.gbif 0.1.0

* First release under the name `wsl.gbif` (October 2022), hosted on EnviDat
  (10.16904/envidat.352) and GitHub (https://github.com/8Ginette8/wsl.gbif).
* Core functions: `wsl_gbif()` (occurrence download with synonym support),
  `wsl_obs_filter()` (grid-based occurrence thinning), `wsl_taxonomy()`
  (GBIF backbone taxonomy lookup), `wsl_doi()` (GBIF-derived DOI generation),
  `make_tiles()` (geographic tiling for `rgbif`).
* Credential-free GBIF download with 13 post-processing filters via
  `CoordinateCleaner`.
* Dynamic moving-window tiling for datasets exceeding 100,000 observations.
