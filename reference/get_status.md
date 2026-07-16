# Retrieve Taxonomy and IUCN Status from GBIF

Query the GBIF backbone taxonomy for a species name, report the matched
accepted name and synonyms, and return the associated IUCN status when
available.

## Usage

``` r
get_status(
  sp_name = NULL,
  search = TRUE,
  rank = NULL,
  phylum = NULL,
  class = NULL,
  order = NULL,
  family = NULL,
  conf_match = 80,
  level = c("accepted", "children", "all"),
  verbose = TRUE
)
```

## Arguments

- sp_name:

  Character. String with the species name. Scientific names at
  genus-species level are expected; fuzzy matching is available when
  `search = FALSE`.

- search:

  Logical. If `TRUE` (default), use a strict GBIF backbone search and
  keep only species-, subspecies-, or variety-level matches. If `FALSE`,
  use a more permissive search and optionally rely on `rank`, `phylum`,
  `class`, `order`, and `family` to resolve ambiguous matches.

- rank:

  Character. String giving the preferred rank to keep: `"SPECIES"`,
  `"SUBSPECIES"`, or `"VARIETY"`. When `NULL`, rank priority is inferred
  from `sp_name`.

- phylum:

  Optional character. Phylum used to disambiguate alternative GBIF
  matches. Particularly useful for hemihomonyms.

- class:

  Optional character. Class used to disambiguate alternative GBIF
  matches.

- order:

  Optional character. Order used to disambiguate alternative GBIF
  matches.

- family:

  Optional character. Family used to disambiguate alternative GBIF
  matches.

- conf_match:

  Numeric. Confidence threshold between 0 and 100 for the GBIF backbone
  match. Default is `80`.

- level:

  Character. Controls how much of the GBIF taxon concept is returned.
  `"accepted"` (default) returns the accepted name and its synonyms.
  `"children"` additionally includes infra-specific taxa (subspecies,
  varieties) as `CHILDREN` rows, mirroring the scope of
  [`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md).
  `"all"` further adds alternative name string representations
  (orthographic variants, etc.) as `RELATED` rows, which are not used by
  [`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
  and are intended for taxonomic inspection only.

- verbose:

  Logical. Should status messages be printed to the console when no
  match is found? Default is `TRUE`.

## Value

A data frame with the columns `canonicalName`, `rank`, `gbif_key`,
`scientificName`, `gbif_status`, `Genus`, `Family`, `Order`, `Class`,
`Phylum`, `IUCN_status`, and `sp_nameMatch`. The row with
`gbif_status = "ACCEPTED"` identifies the accepted GBIF taxon concept
used by
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md),
while `sp_nameMatch` marks the row that most closely matches the
submitted input name. Additional rows carry `gbif_status = "CHILDREN"`
(infra-specific taxa, when `level = "children"` or `level = "all"`) or
`gbif_status = "RELATED"` (alternative name representations, when
`level = "all"`).

## Details

Use `get_status()` before
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
when you want to inspect how GBIF resolves an input name, which name is
treated as the accepted backbone taxon, and which synonyms are included
in the taxon concept used for occurrence retrieval.

When `level = "accepted"`, the returned rows correspond to the accepted
name and synonyms linked to the accepted GBIF taxon key. When
`level = "children"`, infra-specific taxa that
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
also retrieves are added, giving a complete picture of the occurrence
download scope. When `level = "all"`, alternative name representations
are further included for taxonomic review.

IUCN Red List status is retrieved from GBIF for the accepted taxon key
only. GBIF does not currently store subspecies-level IUCN assessments,
so the `IUCN_status` column is uniform across all rows. Note that some
subspecies have independent IUCN assessments (e.g. *Panthera tigris
sumatrae* is Critically Endangered) that are not reflected here; consult
the IUCN Red List directly for infra-specific status.

## References

Chamberlain, S., Oldoni, D., & Waller, J. (2022). rgbif: interface to
the global biodiversity information facility API.
[doi:10.5281/zenodo.6023735](https://doi.org/10.5281/zenodo.6023735)

## See also

[`get_gbif`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)()
to download occurrences for the accepted taxon concept returned here;
the rgbif package for more general approaches to querying the GBIF
backbone taxonomy.

## Examples

``` r
# \donttest{
# --- 1. Default: accepted name + synonyms only ---
tax <- get_status("Panthera tigris")
tax
#>           canonicalName    rank gbif_key                   scientificName
#> 5219416 Panthera tigris SPECIES  5219416 Panthera tigris (Linnaeus, 1758)
#> 4969819    Felis tigris SPECIES  4969819      Felis tigris Linnaeus, 1758
#>         gbif_status    Genus  Family     Order    Class   Phylum IUCN_status
#> 5219416    ACCEPTED Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#> 4969819     SYNONYM Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#>         sp_nameMatch
#> 5219416        INPUT
#> 4969819        EXACT

# --- 2. Include infra-specific taxa (subspecies, varieties) ---
# These are the same keys retrieved by get_gbif()
tax_ch <- get_status("Panthera tigris", level = "children")
tax_ch
#>                     canonicalName       rank gbif_key
#> 5219416           Panthera tigris    SPECIES  5219416
#> 4969819              Felis tigris    SPECIES  4969819
#> 5219420   Panthera tigris altaica SUBSPECIES  5219420
#> 5219425 Panthera tigris amoyensis SUBSPECIES  5219425
#> 5219419    Panthera tigris balica SUBSPECIES  5219419
#> 5219424  Panthera tigris corbetti SUBSPECIES  5219424
#> 5219422  Panthera tigris sondaica SUBSPECIES  5219422
#> 5219418  Panthera tigris sumatrae SUBSPECIES  5219418
#> 7059276    Panthera tigris tigris SUBSPECIES  7059276
#> 5219423   Panthera tigris virgata SUBSPECIES  5219423
#>                                       scientificName gbif_status    Genus
#> 5219416             Panthera tigris (Linnaeus, 1758)    ACCEPTED Panthera
#> 4969819                  Felis tigris Linnaeus, 1758     SYNONYM Panthera
#> 5219420       Panthera tigris altaica Temminck, 1844    CHILDREN Panthera
#> 5219425 Panthera tigris amoyensis (Hilzheimer, 1905)    CHILDREN Panthera
#> 5219419         Panthera tigris balica Schwarz, 1912    CHILDREN Panthera
#> 5219424         Panthera tigris corbetti Mazak, 1968    CHILDREN Panthera
#> 5219422      Panthera tigris sondaica Temminck, 1844    CHILDREN Panthera
#> 5219418        Panthera tigris sumatrae Pocock, 1929    CHILDREN Panthera
#> 7059276                       Panthera tigris tigris    CHILDREN Panthera
#> 5219423      Panthera tigris virgata (Illiger, 1815)    CHILDREN Panthera
#>          Family     Order    Class   Phylum IUCN_status sp_nameMatch
#> 5219416 Felidae Carnivora Mammalia Chordata  ENDANGERED        INPUT
#> 4969819 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219420 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219425 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219419 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219424 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219422 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219418 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 7059276 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5219423 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT

# --- 3. Include alternative name string representations (inspection only) ---
# RELATED rows are not used by get_gbif()
tax_rel <- get_status("Panthera tigris", level = "all")
tax_rel
#>                canonicalName       rank gbif_key
#> 1            Panthera tigris    SPECIES  5219416
#> 2               Felis tigris    SPECIES  4969819
#> 3    Panthera tigris altaica SUBSPECIES  5219420
#> 4  Panthera tigris amoyensis SUBSPECIES  5219425
#> 5     Panthera tigris balica SUBSPECIES  5219419
#> 6   Panthera tigris corbetti SUBSPECIES  5219424
#> 7   Panthera tigris sondaica SUBSPECIES  5219422
#> 8   Panthera tigris sumatrae SUBSPECIES  5219418
#> 9     Panthera tigris tigris SUBSPECIES  7059276
#> 10   Panthera tigris virgata SUBSPECIES  5219423
#> 11           Panthera tigris    SPECIES  5219416
#> 13           Panthera tigris    SPECIES  5219416
#> 14           Panthera tigris    SPECIES  5219416
#> 15           Panthera tigris    SPECIES  5219416
#> 17              Felis tigris    SPECIES  4969819
#> 18              Felis tigris    SPECIES  4969819
#> 19              Felis tigris    SPECIES  4969819
#> 20   Panthera tigris altaica SUBSPECIES  5219420
#> 21   Panthera tigris altaica SUBSPECIES  5219420
#> 22  Panthera tigris altaicus SUBSPECIES  5219420
#> 23   Panthera tigris altaica SUBSPECIES  5219420
#> 24   Panthera tigris altaica SUBSPECIES  5219420
#> 26 Panthera tigris amoyensis SUBSPECIES  5219425
#> 27 Panthera tigris amoyensis SUBSPECIES  5219425
#> 28 Panthera tigris amoyensis SUBSPECIES  5219425
#> 29 Panthera tigris amoyensis SUBSPECIES  5219425
#> 30 Panthera tigris amoyensis SUBSPECIES  5219425
#> 32    Panthera tigris balica SUBSPECIES  5219419
#> 34    Panthera tigris balica SUBSPECIES  5219419
#> 35  Panthera tigris corbetti SUBSPECIES  5219424
#> 36  Panthera tigris corbetti SUBSPECIES  5219424
#> 37  Panthera tigris corbetti SUBSPECIES  5219424
#> 39  Panthera tigris sondaica SUBSPECIES  5219422
#> 40  Panthera tigris sondaica SUBSPECIES  5219422
#> 42  Panthera tigris sondaica SUBSPECIES  5219422
#> 43  Panthera tigris sondaica SUBSPECIES  5219422
#> 44  Panthera tigris sumatrae SUBSPECIES  5219418
#> 45  Panthera tigris sumatrae SUBSPECIES  5219418
#> 47  Panthera tigris sumatrae SUBSPECIES  5219418
#> 48    Panthera tigris tigris SUBSPECIES  7059276
#> 49    Panthera tigris tigris SUBSPECIES  7059276
#> 51    Panthera tigris tigris SUBSPECIES  7059276
#> 52   Panthera tigris virgata SUBSPECIES  5219423
#> 53   Panthera tigris virgata SUBSPECIES  5219423
#> 54   Panthera tigris virgata SUBSPECIES  5219423
#> 55   Panthera tigris virgata SUBSPECIES  5219423
#> 56   Panthera tigris virgata SUBSPECIES  5219423
#>                                         scientificName gbif_status    Genus
#> 1                     Panthera tigris (Linnaeus, 1758)    ACCEPTED Panthera
#> 2                          Felis tigris Linnaeus, 1758     SYNONYM Panthera
#> 3               Panthera tigris altaica Temminck, 1844    CHILDREN Panthera
#> 4         Panthera tigris amoyensis (Hilzheimer, 1905)    CHILDREN Panthera
#> 5                 Panthera tigris balica Schwarz, 1912    CHILDREN Panthera
#> 6                 Panthera tigris corbetti Mazak, 1968    CHILDREN Panthera
#> 7              Panthera tigris sondaica Temminck, 1844    CHILDREN Panthera
#> 8                Panthera tigris sumatrae Pocock, 1929    CHILDREN Panthera
#> 9                               Panthera tigris tigris    CHILDREN Panthera
#> 10             Panthera tigris virgata (Illiger, 1815)    CHILDREN Panthera
#> 11                     Panthera tigris (Linnaeus 1758)     RELATED Panthera
#> 13                    Panthera tigris (‎Linnaeus‎, ‎1758‎)     RELATED Panthera
#> 14                                     Panthera tigris     RELATED Panthera
#> 15                      Panthera tigris Linnaeus, 1758     RELATED Panthera
#> 17                                        Felis tigris     RELATED Panthera
#> 18                          Felis tigris Linnaeus 1758     RELATED Panthera
#> 19                       Felis tigris (Linnaeus, 1758)     RELATED Panthera
#> 20       Panthera tigris subsp. altaica Temminck, 1844     RELATED Panthera
#> 21            Panthera tigris altaica (Temminck, 1844)     RELATED Panthera
#> 22             Panthera tigris altaicus Temminck, 1844     RELATED Panthera
#> 23        Panthera tigris subsp. altaica Temminck 1844     RELATED Panthera
#> 24                             Panthera tigris altaica     RELATED Panthera
#> 26 Panthera tigris subsp. amoyensis (Hilzheimer, 1905)     RELATED Panthera
#> 27          Panthera tigris amoyensis Hilzheimer, 1905     RELATED Panthera
#> 28    Panthera tigris subsp. amoyensis Hilzheimer 1905     RELATED Panthera
#> 29                           Panthera tigris amoyensis     RELATED Panthera
#> 30   Panthera tigris subsp. amoyensis Hilzheimer, 1905     RELATED Panthera
#> 32         Panthera tigris subsp. balica Schwarz, 1912     RELATED Panthera
#> 34               Panthera tigris balica Schwarz, 1912.     RELATED Panthera
#> 35         Panthera tigris subsp. corbetti Mazák, 1968     RELATED Panthera
#> 36                            Panthera tigris corbetti     RELATED Panthera
#> 37         Panthera tigris subsp. corbetti Mazak, 1968     RELATED Panthera
#> 39      Panthera tigris subsp. sondaica Temminck, 1844     RELATED Panthera
#> 40                            Panthera tigris sondaica     RELATED Panthera
#> 42            Panthera tigris sondaica Temminck, 1844.     RELATED Panthera
#> 43    Panthera tigris subsp. sondaica (Temminck, 1844)     RELATED Panthera
#> 44                            Panthera tigris sumatrae     RELATED Panthera
#> 45        Panthera tigris subsp. sumatrae Pocock, 1929     RELATED Panthera
#> 47                Panthera tigris sumatrae Pocock 1929     RELATED Panthera
#> 48                       Panthera tigris subsp. tigris     RELATED Panthera
#> 49         Panthera tigris subsp. tigris Linnaeus 1758     RELATED Panthera
#> 51             Panthera tigris tigris (Linnaeus, 1758)     RELATED Panthera
#> 52      Panthera tigris subsp. virgata (Illiger, 1815)     RELATED Panthera
#> 53               Panthera tigris virgata Illiger, 1815     RELATED Panthera
#> 54         Panthera tigris subsp. virgata Illiger 1815     RELATED Panthera
#> 55        Panthera tigris subsp. virgata Illiger, 1815     RELATED Panthera
#> 56                             Panthera tigris virgata     RELATED Panthera
#>     Family     Order    Class   Phylum IUCN_status sp_nameMatch
#> 1  Felidae Carnivora Mammalia Chordata  ENDANGERED        INPUT
#> 2  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 3  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 4  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 5  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 6  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 7  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 8  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 9  Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 10 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 11 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 13 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 14 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 15 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 17 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 18 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 19 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 20 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 21 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 22 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 23 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 24 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 26 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 27 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 28 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 29 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 30 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 32 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 34 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 35 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 36 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 37 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 39 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 40 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 42 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 43 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 44 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 45 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 47 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 48 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 49 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 51 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 52 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 53 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 54 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 55 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT
#> 56 Felidae Carnivora Mammalia Chordata  ENDANGERED        EXACT

# --- 4. Fuzzy matching for uncertain names ---
# sp_nameMatch = "VARIANT" for close but non-identical fuzzy matches
get_status("Panthera tigri", search = FALSE)
#>           canonicalName    rank gbif_key                   scientificName
#> 5219416 Panthera tigris SPECIES  5219416 Panthera tigris (Linnaeus, 1758)
#> 4969819    Felis tigris SPECIES  4969819      Felis tigris Linnaeus, 1758
#>         gbif_status    Genus  Family     Order    Class   Phylum IUCN_status
#> 5219416    ACCEPTED Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#> 4969819     SYNONYM Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#>         sp_nameMatch
#> 5219416        INPUT
#> 4969819      VARIANT

# --- 5. Cross-check get_status() keys against get_gbif() output ---
occ <- get_gbif("Panthera tigris", has_xy = TRUE, verbose = FALSE)
valid_keys    <- tax_ch$gbif_key[tax_ch$gbif_status %in% c("ACCEPTED", "CHILDREN")]
returned_keys <- unique(occ$acceptedTaxonKey)

# Keys in get_status() — backbone taxa (ACCEPTED + CHILDREN)
valid_keys
#> [1] "5219416" "5219420" "5219425" "5219419" "5219424" "5219422" "5219418"
#> [8] "7059276" "5219423"

# Keys returned by get_gbif()
returned_keys
#> [1]  7059276  5219416  5219420  5219422  5219418 10568627 10604898

# Are all occurrence keys known backbone keys?
# FALSE indicates non-backbone entries (e.g. BOLD sequences) are present
all(returned_keys %in% valid_keys)
#> [1] FALSE

# Inspect non-backbone records — typically MATERIAL_SAMPLE, excluded by default
extra_keys <- returned_keys[!returned_keys %in% valid_keys]
extra_keys
#> [1] 10568627 10604898
occ[occ$acceptedTaxonKey %in% extra_keys,
    c("acceptedTaxonKey", "scientificName", "basisOfRecord")]
#>      acceptedTaxonKey scientificName   basisOfRecord
#> 4307         10568627   BOLD:AAD6820 MATERIAL_SAMPLE
#> 4557         10568627   BOLD:AAD6820 MATERIAL_SAMPLE
#> 4573         10568627   BOLD:AAD6820 MATERIAL_SAMPLE
#> 5292         10604898   BOLD:AAC3048 MATERIAL_SAMPLE
#> 5302         10604898   BOLD:AAC3048 MATERIAL_SAMPLE
#> 5303         10604898   BOLD:AAC3048 MATERIAL_SAMPLE

# Note: not all valid_keys need to appear in returned_keys — some subspecies
# may have no records in GBIF at all (e.g. extinct taxa) or may have been
# removed by get_gbif() filtering options (e.g. has_xy, basis, grain)

# --- 6. Input name flagged correctly regardless of how GBIF resolves it ---
# sp_nameMatch = "INPUT" marks the row closest to the submitted name
tax_ch[tax_ch$sp_nameMatch == "INPUT", ]
#>           canonicalName    rank gbif_key                   scientificName
#> 5219416 Panthera tigris SPECIES  5219416 Panthera tigris (Linnaeus, 1758)
#>         gbif_status    Genus  Family     Order    Class   Phylum IUCN_status
#> 5219416    ACCEPTED Panthera Felidae Carnivora Mammalia Chordata  ENDANGERED
#>         sp_nameMatch
#> 5219416        INPUT

# }
```
