# Package index

## GBIF retrieval

Download and inspect GBIF occurrences

- [`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
  : Download and Filter GBIF Occurrences for Spatial Analyses
- [`get_gbif_count()`](https://8ginette8.github.io/gbif.range/reference/get_gbif_count.md)
  : Count GBIF Occurrences Before Downloading Data
- [`get_status()`](https://8ginette8.github.io/gbif.range/reference/get_status.md)
  : Retrieve Taxonomy and IUCN Status from GBIF
- [`get_doi()`](https://8ginette8.github.io/gbif.range/reference/get_doi.md)
  : Create a DOI for GBIF-derived datasets
- [`make_tiles()`](https://8ginette8.github.io/gbif.range/reference/make_tiles.md)
  : Create Tiled GBIF Geometry Queries

## Range inference

Build and evaluate ecoregion-constrained range maps

- [`get_range()`](https://8ginette8.github.io/gbif.range/reference/get_range.md)
  : Create Species Range Maps from Occurrences and Ecoregions
- [`cv_range()`](https://8ginette8.github.io/gbif.range/reference/cv_range.md)
  : Evaluate a Range Map by Cross-Validation
- [`evaluate_range()`](https://8ginette8.github.io/gbif.range/reference/evaluate_range.md)
  : Evaluate Range Maps Against Independent Validation Data
- [`make_blocks()`](https://8ginette8.github.io/gbif.range/reference/make_blocks.md)
  : Split Data into Approximately Balanced Folds

## Ecoregions

Load, download, or create ecoregion layers

- [`read_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)
  : Read an Ecoregion Layer
- [`make_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/make_ecoreg.md)
  : Build Custom Ecoregions from Environmental Layers
- [`get_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/get_ecoreg.md)
  : Download Ecoregion Layers
- [`check_and_get_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/check_and_get_ecoreg.md)
  : Check for a Local Ecoregion Layer and Download It if Needed
- [`ecoreg_list`](https://8ginette8.github.io/gbif.range/reference/ecoreg_list.md)
  : Metadata for Downloadable Ecoregion Layers

## Occurrence filtering

Thin or aggregate occurrence records

- [`obs_filter()`](https://8ginette8.github.io/gbif.range/reference/obs_filter.md)
  : Filter GBIF Records by Grid Cell

## Large-scale batch workflow

Disk-based processing for large GBIF exports

- [`split_gbif_by_species()`](https://8ginette8.github.io/gbif.range/reference/split_gbif_by_species.md)
  : Split a Downloaded GBIF Table into One File per Species

- [`species_csvs_to_ranges()`](https://8ginette8.github.io/gbif.range/reference/species_csvs_to_ranges.md)
  : Build Range Maps from Per-Species GBIF Files Saved on Disk

- [`read_range_rds()`](https://8ginette8.github.io/gbif.range/reference/read_range_rds.md)
  :

  Read a Range File Saved by
  [`species_csvs_to_ranges()`](https://8ginette8.github.io/gbif.range/reference/species_csvs_to_ranges.md)

## Data

Datasets bundled with the package

- [`area_data`](https://8ginette8.github.io/gbif.range/reference/area_data.md)
  : Range Area Comparison Data

## Package

- [`gbif.range`](https://8ginette8.github.io/gbif.range/reference/gbif.range.md)
  [`gbif.range-package`](https://8ginette8.github.io/gbif.range/reference/gbif.range.md)
  [`package-gbif.range`](https://8ginette8.github.io/gbif.range/reference/gbif.range.md)
  : gbif.range: Tools for GBIF Retrieval and Ecoregion-Based Species
  Range Mapping
