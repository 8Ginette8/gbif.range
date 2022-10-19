---
title: 'wsl.gbif - A R package to extract and compile GBIF datasets for large-scale analyses'
tags:
  - R
  - Ecology
  - Geographic Information Systems
  - Species records
  - Scientific names
  - Global Biodiversity Information Facility
authors:
  - name: Yohann Chauvier
    orcid: 0000-0001-9399-3192
    corresponding: true
    affiliation: 1 
  - name: Patrice Descombes
    orcid: 0000-0002-3760-9907
    affiliation: "2, 3"
  - name: Michael P. Nobis
    affiliation: 4
  - name:  Katalin Csilléry
    orcid: 0000-0003-0039-9296
    affiliation: 1
affiliations:
 - name: Evolutionary Genetics, Biodiversity and Conservation Biology, Swiss Federal Research Institute (WSL), 8903 Birmensdorf, Switzerland
   index: 1
 - name: Musée et Jardins botaniques cantonaux, Av. de Cour 14B, 1007 Lausanne, Switzerland
   index: 2
 - name: École Polytechnique Fédérale de Lausanne, Rte Cantonale, 1015 Lausanne, Switzerland
   index: 3
 - name: Dynamic Macroecology, Land Change Science, Swiss Federal Research Institute (WSL), 8903 Birmensdorf, Switzerland
   index: 4
date: 18.10.2022
bibliography: paper.bib

---

# Summary

The Global Biodiversity Information facility (GBIF) is a large public data repository
inventorying georeferenced species observations from all taxonomic groups, with an
increasing number of public, associative and private contributors that updates this
database everyday by uploading their own datasets. The GBIF repository is freely
accessible (https://www.gbif.org/), includes > two billions observations (data from
2022) and is frequently used by a vast community of ecologists, modellers, researchers
and practitioners to adress ecological, geological and environmental questions
whose answers are becoming increasingly crucial in our current global change contetx.
Since rencently, retrieving GBIF species information may easily be done
via new packages and libraries, e.g. *rgbif* on R (Chamberlain & al. 2022) or
*pygbif* on python (Blumentrath & Kudrnovsky 2016). However, the given accessibility
of these libraries are impended by their technicalities, various functions and
parameters, and limitations in the number of species observations a user may download
at once (< 100,000 records). Here we present *wsl.gbif*, a R package which simplifies
the extraction and compilation of large GBIF datasets by providing a user-friendly main
search function *wsl_gbif* that automatically retrieve species records linked to several
scientific names, bypasses the download hard limit of 100,000 observations, and filter
observations quality via several simple parameters.

# Statement of need

While the *rgbif* R package offers great solutions to efficiently extract observational
datasets from GBIF, the package function is not meant to automatically download large GBIF
datasets of millions of records. The downlaod hard limit is set to 100,000 observations which
provide complications for large-scale studies or analyses that focus on several regions of the
world and many taxa. Also, *rgbif* includes several search functions including many parameters
that impend everyday users to easily extract the same number of species records found on the main
website (https://www.gbif.org/). wsl.gbif therefore aims at easing the workflow of retrieving,
compilling GBIF observations at large spatial scale for all species accepted names and synonyms,
and selecting them according to the specific scale of a spatial analysis. The package has and is
currently used in several publications and ongoing project (Chauvier & al. 2021, 2022; ...)

On the one hand, *wsl.gbif* allows the whole observations of a given species name (accepted and
synonym names) to be retrieved, and improves the data accessibility of the *rgbif* R package
(https://cran.r-project.org/web/packages/rgbif/index.html). The user download hard limit of
rgbif is a maximum of 100,000 of species observations in one go (https://www.gbif.org/developer/occurrence,
'Searching' header). This impends the accessibility to the GBIF database when large observational
datasets for many species and regions of the world are needed, specifically in scientifc fields
related to macrocology, modelling or satelite imagery. wsl.gbif therefore bypasses this limit
using geographic parameters from the *rgbif* and *terra* R packages and by adopting a dynamic moving
window process allowing the user's study area of interest to be automatically fragmented in several
tiles that always include < 100,000 observations.

On the other hand, *wsl.gbif* implements easy to use preliminary filtering options implemented during
the download so that the user saves some post-processing time in data cleaning. The filetring
otpions (n = 11) may be chosen independently, and two of them are based on the *CoordinateCleaner* R
package (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html). It is important to
note that, although a strong records filetring may be undertaken with wsl.gbif, *CoordinateCleaner*
includes a larger variety of options that should be checked and applied on *wsl.gbif* outputs.

Finally, *wsl.gbif* offers a set of additional very useful functions meant to be used for large-scale
studies using GBIF observations. (1) *wsl_taXnames* generates, based on a given species name, a list
of all its scientific names (accepted, synonyms, children and related) found in the GBIF backbone
taxonomy. The function allows therefore taxonomy correspondency to be made between different species
and sub-species to potentially merge their records, but also permits efficient ways of linking external
data of a species which is named differently across databases. (2) Whereas the 'grain' parameter
in *wsl_gbif* allows GBIF observations to be filtered according to a certain spatial precision,
*wsl_obs_filter* accepts as input a *wsl_gbif* output (one or several species) and filter the observations
according to a specific given grid resolution (one observation per pixel grid kept). (3) wsl_tiles is
a function that may be used to generate a set of n geometry arguments POLYGON() based on a given
geographic extent. This function is meant to help users who want to use the *rgbif* R package and its
parameter *geometry* that uses a POLYGON() argument.

# Examples

Load the library

``` r
library(wsl.gbif)
```

Let's download worldwide the observations of Panthera tigris:

``` r
obs.pt = wsl_gbif("Panthera tigris")
```

Or simply get all its scientific names (accepted and synonyms) from the GBIF backbone taxonomy:

``` r
wsl_taXnames("Panthera tigris",all=FALSE)
```

Same may be done with Delphinus delphis (a species with > 100'00 observations)

``` r
obs. dd = wsl_gbif("Delphinus delphis")
wsl_taXnames("Delphinus delphis",all=TRUE) # Here the list is longer because 'all=TRUE' includes every names (even doubtful)
```


# Acknowledgements

...

# References

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity information facility API. <a href="https://doi.org/10.5281/zenodo.6023735">10.5281/zenodo.6023735</a>

Blumentrath, S., & Kudrnovsky, H. (2016). v. in. pygbif-Search and import GBIF species distribution data.

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>

Chauvier, Y., Descombes, P., Guéguen, M., Boulangeat, L., Thuiller, W., & Zimmermann, N. E. (2022). Resolution in species distribution models shapes spatial patterns of plant multifaceted diversity. Ecography, 2022 (10), e05973. <a href="https://doi.org/10.1111/ecog.05973">10.1111/ecog.05973</a>

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5), 744-751. <a href="https://doi.org/10.1111/2041-210X.13152">10.1111/2041-210X.13152</a>

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). <a href="https://cran.r-project.org/web/packages/terra/index.html">Terra - CRAN</a>



