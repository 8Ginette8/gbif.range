# Range Area Comparison Data

A dataset containing range-area estimates for species, comparing
`gbif.range`-derived polygons with corresponding IUCN ranges.

## Format

A data frame with 3 columns:

- species:

  Character string with the scientific name of the species.

- gbif.range.km2:

  Numeric area in square kilometers calculated from `gbif.range`
  polygons.

- iucn.range.km2:

  Numeric area in square kilometers calculated from IUCN range data.

## Source

Calculated with the `gbif.range` package.
