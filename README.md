# wsl.gbif

wsl.gbif aims at easing the workflow of retrieving GBIF observations at large spatial scale for all species accepted names and synonyms, and filtering them according to the specific scale of a spatial analysis.

On the one hand, wsl.gbif allows the whole observations of a given species name (accepted and synonym names) to be retrieved, and improves the data accessibility of the rgbif R package (https://cran.r-project.org/web/packages/rgbif/index.html). The user download limit of rgbif is a maximum of 100'000 of species observations in one go. This impends the accessibility to the GBIF database when large observational datasets for many species and regions of the world are neededed. wsl.gbif therefore bypasses this limit using geographic parameters from the rgbig R package and by adopting a dynamic moving window process allowing the user's study area of interest to be fragmented in several tiles that always include < 100'000 observations.

On the other hand, wsl.gbif implements easy to use preliminary filtering options implemented during the download so that the user saves some post-processing time in data cleaning. The filetring otpions (n = 11) may be chosen independently, and two of them are based on the CoordinateCleaner R package (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html). Additionally, wsl.gbif offers a second function that generates, based on a given species name, a list of all its scientific names (accepted, synonyms, children and related) found in the GBIF backbone taxonomy.

## Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("8Ginette8/wsl.gbif")
```
The package will soon be the object of a peer review process, and will therefore be published on CRAN.

## Example

Let's download worldwide the observations of Panthera tigris:

``` r
wsl_gbif("Panthera tigris")
```

Or simply get all its scientific names (accepted and synonyms) from the GBIF backbone taxonomy:

``` r
wsl_taXnames("Cypripedium calceolus",all=FALSE)
```

Same may be done with Delphinus delphis (a species with > 100'00 observations)

``` r
wsl_gbif("Delphinus delphis")
wsl_taXnames("Delphinus delphis",all=TRUE) # Here the list is long because all=TRUE inlcudes every names (even doubtful)
```


## References

@article{chauvier2021influence,
  title={Influence of climate, soil, and land cover on plant species distribution in the European Alps},
  author={Chauvier, Yohann and Thuiller, Wilfried and Brun, Philipp and Lavergne, S{\'e}bastien and Descombes, Patrice and Karger, Dirk N and Renaud, Julien and Zimmermann, Niklaus E},
  journal={Ecological monographs},
  volume={91},
  number={2},
  pages={e01433},
  year={2021},
  publisher={Wiley Online Library}
}
