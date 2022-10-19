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
via new packages and libraries, e.g. rgbif on R (Chamberlain & al. 2022) or
pygbif on python (Blumentrath & Kudrnovsky 2016). However, the given accessibility
of these libraries are impended by their technicalities, various functions and
parameters, and limitations in the number of species observations a user may download
at once (< 100,000 records). Here we present wsl.gbif, a R package which simplifies
the extraction and compilation of large GBIF datasets by providing a user-friendly main
search function (wsl_gbif) that automatically retrieve species records linked to several scientific
names, bypasses the download hard limit of 100,000 observations, and filter observations
quality via several simple parameters.

# Statement of need

While the rgbif package offers great solutions to efficiently extract observational datasets
from GBIF, the package function is not meant to automatically download large GBIF datasets
of millions of records. The downlaod hard limit is set to 100,000 observations which provide
complications for large-scale studies or analyses that focus on several regions of the world
and many taxa. Also, rgbif includes several search functions including many parameters that
impend everyday users to easily extract the same number of species records found on the main
website (https://www.gbif.org/). wsl.gbif therefore aims at easing the workflow of retrieving,
compilling GBIF observations at large spatial scale for all species accepted names and synonyms,
and select them according to the specific scale of a spatial analysis.

On the one hand, wsl.gbif allows the whole observations of a given species name (accepted and
synonym names) to be retrieved, and improves the data accessibility of the rgbif R package
(https://cran.r-project.org/web/packages/rgbif/index.html). The user download hard limit of
rgbif is a maximum of 100,000 of species observations in one go (https://www.gbif.org/developer/occurrence,
'Searching' header). This impends the accessibility to the GBIF database when large observational
datasets for many species and regions of the world are needed, specifically in scientifc fields
related to macrocology, modelling or satelite imagery. wsl.gbif therefore bypasses this limit
using geographic parameters from the rgbif and terra R packages and by adopting a dynamic moving
window process allowing the user's study area of interest to be automatically fragmented in several
tiles that always include < 100,000 observations.

On the other hand, wsl.gbif implements easy to use preliminary filtering options implemented during
the download so that the user saves some post-processing time in data cleaning. The filetring
otpions (n = 11) may be chosen independently, and two of them are based on the CoordinateCleaner R
package (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html). It is important to
note that, although a strong records filetring may be undertaken with wsl.gbif, CoordinateCleaner
includes a larger variety of options that should be checked and applied on wsl.gbif outputs.

Finally, wsl.gbif offers a set of additional very useful functions meant to be used for large-scale
studies. (1) wsl_taXnames generates, based on a given species name, a list of all its scientific names
(accepted, synonyms, children and related) found in the GBIF backbone taxonomy. (2) Whereas the 'grain'
parameter in wsl_gbif allows unprecise GBIF observations to be filtered, wsl_obs_filter accepts as input
a wsl_gbif output (one or several species) and keeps for each species only one observation per pixel of a
defined grid. (3) wsl_tiles is a function that may be used to generate a set of n geometry arguments
POLYGON() based on a given geographic extent. This function is meant to help users who want to use the
rgbif package and its parameter 'geometry' that use a POLYGON() argument.


# Examples

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
