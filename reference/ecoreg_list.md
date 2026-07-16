# Metadata for Downloadable Ecoregion Layers

A named list describing the ecoregion layers that can be downloaded with
[`get_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/get_ecoreg.md).
Each element is itself a list with three character fields: `filename`
(the local save name), `link` (the download URL), and `description` (a
short label).

## Usage

``` r
ecoreg_list
```

## Format

A list with each element containing:

- filename:

  The name to save the downloaded file as.

- link:

  The URL to download the file.

- description:

  A short description of the file.

## Value

A named list of length 4, one element per available ecoregion dataset,
each with `filename`, `link`, and `description` fields as described in
`Format`.

## See also

[`get_ecoreg`](https://8ginette8.github.io/gbif.range/reference/get_ecoreg.md)()
to download the datasets listed here, and
[`read_ecoreg`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)()
to load one directly.

## Examples

``` r
ecoreg_list
#> [[1]]
#> [[1]]$filename
#> [1] "eco_terra"
#> 
#> [[1]]$link
#> [1] "https://www.arcgis.com/sharing/rest/content/items/b1636d640ede4d6ca8f5e369f2dc368b/data"
#> 
#> [[1]]$description
#> [1] "Terrestrial Ecoregions of the World"
#> 
#> 
#> [[2]]
#> [[2]]$filename
#> [1] "eco_fresh"
#> 
#> [[2]]$link
#> [1] "https://www.arcgis.com/sharing/rest/content/items/41aa8662254f43699e792fa194d7b3bb/data"
#> 
#> [[2]]$description
#> [1] "Freshwater Ecoregions of the World"
#> 
#> 
#> [[3]]
#> [[3]]$filename
#> [1] "eco_marine"
#> 
#> [[3]]$link
#> [1] "https://www.arcgis.com/sharing/rest/content/items/903c3ae05b264c00a3b5e58a4561b7e6/data"
#> 
#> [[3]]$description
#> [1] "Marine Ecoregions of the World"
#> 
#> 
#> [[4]]
#> [[4]]$filename
#> [1] "eco_hd_marine"
#> 
#> [[4]]$link
#> [1] "https://datadownload-production.s3.amazonaws.com/WCMC036_MEOW_PPOW_2007_2012_v1.zip"
#> 
#> [[4]]$description
#> [1] "High Definition Marine Ecoregions and Pelagic Provinces of the World (2007, 2012)"
#> 
#> 
```
