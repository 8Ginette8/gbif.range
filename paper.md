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

The Global Biodiversity Information facility (GBIF) is a large public data repository inventorying georeferenced species observations from all taxonomic groups, with an increasing number of public, associative and private contributors that updates this database everyday by uploading their own datasets. The GBIF repository is freely accessible (https://www.gbif.org/), includes > two billions observations (data from 2022) and is frequently used by a vast community of ecologists, modellers, researchers and practitioners to address ecological, geological and environmental questions whose answers are becoming increasingly crucial in our current global change context. Since recently, retrieving GBIF species information may easily be done via new packages and libraries, e.g. *rgbif* on R or *pygbif* on python. However, the given accessibility of these libraries are impended by their technicalities, various functions and parameters, and limitations in the number of species observations a user may download at once when no official data request is made (< 100,000 records). Here we present *wsl.gbif*, a small library which simplifies the extraction and compilation of large GBIF datasets by providing a user-friendly main search function wsl_gbif() that automatically retrieve species records linked to several scientific names, bypasses the download hard limit of 100,000 observations of occ_data() in *rgbif*, and filter observations quality via several simple parameters.

# Statement of need

While the *rgbif* R package offers great solutions to efficiently extract observational datasets from GBIF, the package occ_data() and occ_search() fast functions are not meant to automatically download large GBIF datasets of millions of records. The downlaod hard limit is set to 100,000 observations which provide complications for large-scale studies or analyses that focus on several regions of the world and many taxa. Also, *rgbif* includes several search functions including many parameters that impend everyday users to easily extract the same number of species records found on the main website (https://www.gbif.org/). wsl.gbif therefore aims at easing the workflow of retrieving, compilling GBIF observations at large spatial scale for all species accepted names and synonyms, and selecting them according to the specific scale of a spatial analysis. Our small library has and is currently used in several publications and ongoing project (Chauvier & al. 2021, 2022; ...)

On the one hand, *wsl_gbif()* allows the whole observations of a given species name (accepted and synonym names) to be automatically retrieved, and improves the data accessibility of the rgbif R package (https://cran.r-project.org/web/packages/rgbif/index.html). The user download hard limit of *rgbif* is a maximum of 100,000 of species observations in one go if the easy-to-use interactive functions *occ_search()* and *occ_data()* are used (i.e., if no official download request is made with *occ_download()*, https://www.gbif.org/developer/occurrence). This impends the fast accessibility to the GBIF database when large observational datasets for many species and regions of the world are needed, specifically in scientifc fields related to macrocology, modelling or satelite imagery. *wsl_gbif()* therefore bypasses this limit by intuitively using geographic parameters from the *rgbif* *occ_data()* function and *terra* R packages and adopting a dynamic moving window process allowing the user's study area of interest to be automatically fragmented in several tiles that always include < 100,000 observations.

On the other hand, *wsl_gbif()* implements easy to use preliminary filtering options implemented during the download so that the user saves some post-processing time in data cleaning. 13 filters are available. Two already are set by default in *wsl_gbif()* (*hasCoordinate = TRUE*, *hasGeospatialIssue=FALSE*) and 11 can be chosen independently, including two that are based on the *CoordinateCleaner* R package (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html). It is important to note that, although a strong records filetring may be undertaken with *wsl_gbif()*, *CoordinateCleaner* includes a larger variety of options that should be checked and applied on *wsl.gbif* outputs.

Finally, *wsl.gbif* offers a set of additional very useful functions meant to be used for large-scale studies using GBIF observations:
  - *wsl_taxonomy*: Generates, based on a given species name, a list of all its scientific names
  (accepted, synonyms) found in the GBIF backbone taxonomy to download the data. Children and related
  doubtful names not used to download the data may also be extracted. The function allows therefore taxonomy
  correspondency to be made between different species and sub-species to potentially merge their records,
  but also permits efficient ways of linking external data of a species which is named differently across databases.
  - *wsl_obs_filter*: Whereas the 'grain' parameter in *wsl_gbif()* allows GBIF observations to be filtered
  according to a certain spatial precision, *wsl_obs_filter()* accepts as input a *wsl_gbif()* output (one or
  several species) and filter the observations according to a specific given grid resolution (one observation
  per pixel grid kept). This function allows the user to refine the density of GBIF observations according to
  a defined analysis/study's resolution.
  - *wsl_tiles*: May be used to generate a set of *SpatialExtent* and geometry arguments POLYGON() based on a given
  geographic extent. This function is meant to help users who want to use the *rgbif* R package and its parameter
  *geometry* that uses a POLYGON() argument.
  - *wsl_doi*: A small wrapper of *derived_dataset()* in *rgbif* that simplifies the obtention of a general DOI
  for a set of several gbif species datasets.

# Examples

Let's download worldwide the records of Panthera tigris only based on true observations:

``` r
# Download
wsl.pt = wsl_gbif("Panthera tigris",basis=c("OBSERVATION","HUMAN_OBSERVATION"))

# Plot species records
library(maptools)
data(wrld_simpl)
plot(wrld_simpl)
points(wsl.pt[,c("decimalLongitude","decimalLatitude")],pch=20,col="#238b4550",cex=4)
```

Note that an additional filtering needs here to be done as one observation is found in the US. A lot of tigers are being captive in this country hence the recorded observation. Therefore *CoordinateCleaner* functions should here be considered thereafter.

We can also retrieve all the tiger scientific names (accepted and synonyms) that were used in the download with the GBIF backbone taxonomy. If all = TRUE, additonal children and related doubtful names may also beextracted (not used in *wsl_gbif()*):

``` r
wsl_taxonomy("Panthera tigris",all=FALSE)
```

Same may be done with Delphinus delphis (a species with > 100'00 observations)

``` r
wsl_gbif("Delphinus delphis")
wsl_taXnames("Delphinus delphis",all=TRUE) # Here the list is longer because 'all=TRUE' includes every names (even doubtful)
```

# References

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to the global biodiversity information facility API. <a href="https://doi.org/10.5281/zenodo.6023735">10.5281/zenodo.6023735</a>

Blumentrath, S., & Kudrnovsky, H. (2016). v. in. pygbif-Search and import GBIF species distribution data.

Chauvier, Y., Thuiller, W., Brun, P., Lavergne, S., Descombes, P., Karger, D. N., ... & Zimmermann, N. E. (2021). Influence of climate, soil, and land cover on plant species distribution in the European Alps. Ecological monographs, 91(2), e01433. <a href="https://doi.org/10.1002/ecm.1433">10.1002/ecm.1433</a>

Chauvier, Y., Descombes, P., Guéguen, M., Boulangeat, L., Thuiller, W., & Zimmermann, N. E. (2022). Resolution in species distribution models shapes spatial patterns of plant multifaceted diversity. Ecography, 2022 (10), e05973. <a href="https://doi.org/10.1111/ecog.05973">10.1111/ecog.05973</a>

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., ... & Antonelli, A. (2019). CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5), 744-751. <a href="https://doi.org/10.1111/2041-210X.13152">10.1111/2041-210X.13152</a>

Hijmans, Robert J. "terra: Spatial Data Analysis. R Package Version 1.6-7." (2022). <a href="https://cran.r-project.org/web/packages/terra/index.html">Terra - CRAN</a>



