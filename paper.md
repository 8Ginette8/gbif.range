---
title: 'wsl.gbif - An R package to extract and filter large observation datasets from GBIF'
tags:
  - R
  - Ecology
  - Geographic Information Systems
  - Species records
  - Taxonomy
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

wsl.gbif aims at easing the workflow of retrieving GBIF observations
at large spatial scale for all species accepted names and synonyms,
and filtering them according to the specific scale of a spatial analysis.

On the one hand, wsl.gbif allows the whole observations of a given species
name (accepted and synonym names) to be retrieved, and improves the data
accessibility of the rgbif R package (https://cran.r-project.org/web/packages/rgbif/index.html).
The user download hard limit of rgbif is a maximum of 100,000 of species observations in one go
(https://www.gbif.org/developer/occurrence, 'Searching' header). This impends the accessibility
to the GBIF database when large observational datasets for many species and regions of the world
are needed, specifically in scientifc fields related to macrocology, modelling or satelite imagery.
wsl.gbif therefore bypasses this limit using geographic parameters from the rgbif and terra R
packages and by adopting a dynamic moving window process allowing the user's study area of interest
to be automatically fragmented in several tiles that always include < 100,000 observations.

On the other hand, wsl.gbif implements easy to use preliminary filtering options implemented during
the download so that the user saves some post-processing time in data cleaning. The filetring
otpions (n = 11) may be chosen independently, and two of them are based on the CoordinateCleaner R
package (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html). It is important to
note that, although a strong records filetring may be undertaken with wsl.gbif, CoordinateCleaner
includes a larger variety of options that should be checked and applied on wsl.gbif outputs.

Finally, wsl.gbif offers a second function that generates, based on a given species name, a list of
all its scientific names (accepted, synonyms, children and related) found in the GBIF backbone taxonomy.

# Statement of need

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
