# Read an Ecoregion Layer

Load one of the packaged ecoregion datasets listed in
[`ecoreg_list()`](https://8ginette8.github.io/gbif.range/reference/ecoreg_list.md).
If the requested data are not available locally, they are downloaded
first.

## Usage

``` r
read_ecoreg(ecoreg_name = "eco_terra", save_dir = NULL, format = "SpatVector")
```

## Arguments

- ecoreg_name:

  Character. File name of the ecoregion dataset to load. See
  `ecoreg_list` for available options.

- save_dir:

  Character. irectory where the downloaded files are stored.

- format:

  `"SpatVector"` (default) or `"sf"`.

## Value

An ecoregion layer as a `SpatVector` or `sf` object, depending on
`format`.

## Details

Four datasets are currently available:

\(1\) `eco_terra` for terrestrial species, based on The Nature
Conservancy version adapted from Olson et al. (2001).

\(2\) `eco_marine` for marine species, based on Spalding et al. (2007,
2012).

\(3\) `eco_hd_marine`, a higher-resolution marine version.

\(4\) `eco_fresh` for freshwater species, based on Abell et al. (2008).

## References

Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I., Strand,
H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F., Ricketts, T. H.,
Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao, P., Kassem, K. R.
2001. Terrestrial ecoregions of the world: a new map of life on Earth.
BioScience 51(11):933-938. doi: 10.1641/0006-3568(2001)051

The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
Units. GIS layers developed by The Nature Conservancy with multiple
partners, combined from Olson et al. (2001), Bailey 1995 and Wiken 1986.
Cambridge (UK): The Nature Conservancy.

Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A.
Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al
Lombana, Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer
Molnar, Cheri A. Recchia, James Robertson, Marine Ecoregions of the
World: A Bioregionalization of Coastal and Shelf Areas, BioScience,
Volume 57, Issue 7, July 2007, Pages 573–583. doi: 10.1641/B570707

Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
Pelagic provinces of the world: a biogeographic classification of the
world’s surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
doi: 10.1016/j.ocecoaman.2011.12.016

The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
of the World. GIS layers developed by The Nature Conservancy with
multiple partners, combined from Spalding et al. (2007) and Spalding et
al. (2012). Cambridge (UK): The Nature Conservancy.

Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice
Kottelat, Nina Bogutskaya, Brian Coad, Nick Mandrak, Salvador Contreras
Balderas, William Bussing, Melanie L. J. Stiassny, Paul Skelton, Gerald
R. Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf,
James Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel,
Eric Wikramanayake, David Olson, Hugo L. López, Roberto E. Reis, John G.
Lundberg, Mark H. Sabaj Pérez, Paulo Petry, Freshwater Ecoregions of the
World: A New Map of Biogeographic Units for Freshwater Biodiversity
Conservation, BioScience, Volume 58, Issue 5, May 2008, Pages 403–414.
doi: 10.1641/B580507

## Examples

``` r
if (FALSE) { # \dontrun{
shp_eco_terra <- read_ecoreg("eco_terra")
plot(shp_eco_terra)
} # }
```
