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
tax <- get_status("Ailuropoda melanoleuca")
tax
#>                  canonicalName    rank gbif_key
#> 2433399 Ailuropoda melanoleuca SPECIES  2433399
#> 9387176 Aeluropus melanoleucus SPECIES  9387176
#> 9379787 Ailuropus melanoleucus SPECIES  9379787
#> 7888034     Ursus melanoleucus SPECIES  7888034
#>                               scientificName gbif_status      Genus  Family
#> 2433399 Ailuropoda melanoleuca (David, 1869)    ACCEPTED Ailuropoda Ursidae
#> 9387176 Aeluropus melanoleucus (David, 1869)     SYNONYM Ailuropoda Ursidae
#> 9379787 Ailuropus melanoleucus (David, 1869)     SYNONYM Ailuropoda Ursidae
#> 7888034       Ursus melanoleucus David, 1869     SYNONYM Ailuropoda Ursidae
#>             Order    Class   Phylum IUCN_status sp_nameMatch
#> 2433399 Carnivora Mammalia Chordata  VULNERABLE        INPUT
#> 9387176 Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 9379787 Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 7888034 Carnivora Mammalia Chordata  VULNERABLE        EXACT

# --- 2. Include infra-specific taxa (subspecies, varieties) ---
# These are the same keys retrieved by get_gbif()
tax_ch <- get_status("Ailuropoda melanoleuca", level = "children")
tax_ch
#>                                  canonicalName       rank gbif_key
#> 2433399                 Ailuropoda melanoleuca    SPECIES  2433399
#> 9387176                 Aeluropus melanoleucus    SPECIES  9387176
#> 9379787                 Ailuropus melanoleucus    SPECIES  9379787
#> 7888034                     Ursus melanoleucus    SPECIES  7888034
#> 9087131     Ailuropoda melanoleuca melanoleuca SUBSPECIES  9087131
#> 9460321 Ailuropoda melanoleuca wulingshanensis SUBSPECIES  9460321
#>                                                     scientificName gbif_status
#> 2433399                       Ailuropoda melanoleuca (David, 1869)    ACCEPTED
#> 9387176                       Aeluropus melanoleucus (David, 1869)     SYNONYM
#> 9379787                       Ailuropus melanoleucus (David, 1869)     SYNONYM
#> 7888034                             Ursus melanoleucus David, 1869     SYNONYM
#> 9087131                         Ailuropoda melanoleuca melanoleuca    CHILDREN
#> 9460321 Ailuropoda melanoleuca wulingshanensis (Wang et al., 1982)    CHILDREN
#>              Genus  Family     Order    Class   Phylum IUCN_status sp_nameMatch
#> 2433399 Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE        INPUT
#> 9387176 Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 9379787 Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 7888034 Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 9087131 Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 9460321 Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE        EXACT

# --- 3. Include alternative name string representations (inspection only) ---
# RELATED rows are not used by get_gbif()
tax_rel <- get_status("Ailuropoda melanoleuca", level = "all")
tax_rel
#>                             canonicalName       rank gbif_key
#> 1                  Ailuropoda melanoleuca    SPECIES  2433399
#> 2                  Aeluropus melanoleucus    SPECIES  9387176
#> 3                  Ailuropus melanoleucus    SPECIES  9379787
#> 4                      Ursus melanoleucus    SPECIES  7888034
#> 5      Ailuropoda melanoleuca melanoleuca SUBSPECIES  9087131
#> 6  Ailuropoda melanoleuca wulingshanensis SUBSPECIES  9460321
#> 8                  Ailuropoda melanoleuca    SPECIES  2433399
#> 9                  Ailuropoda melanoleuca    SPECIES  2433399
#> 10                 Ailuropoda melanoleuca    SPECIES  2433399
#> 14                      Ursus melanoleuca    SPECIES  7888034
#> 15                     Ursus melanoleucus    SPECIES  7888034
#> 16                     Ursus melanoleucus    SPECIES  7888034
#> 18 Ailuropoda melanoleuca wulingshanensis SUBSPECIES  9460321
#>                                                       scientificName
#> 1                               Ailuropoda melanoleuca (David, 1869)
#> 2                               Aeluropus melanoleucus (David, 1869)
#> 3                               Ailuropus melanoleucus (David, 1869)
#> 4                                     Ursus melanoleucus David, 1869
#> 5                                 Ailuropoda melanoleuca melanoleuca
#> 6         Ailuropoda melanoleuca wulingshanensis (Wang et al., 1982)
#> 8                                Ailuropoda melanoleuca (David 1869)
#> 9                                 Ailuropoda melanoleuca David, 1869
#> 10                                            Ailuropoda melanoleuca
#> 14                                                 Ursus melanoleuca
#> 15                                     Ursus melanoleucus David 1869
#> 16                                  Ursus melanoleucus (David, 1869)
#> 18 Ailuropoda melanoleuca subsp. wulingshanensis (Wang et al., 1982)
#>    gbif_status      Genus  Family     Order    Class   Phylum IUCN_status
#> 1     ACCEPTED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 2      SYNONYM Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 3      SYNONYM Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 4      SYNONYM Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 5     CHILDREN Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 6     CHILDREN Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 8      RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 9      RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 10     RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 14     RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 15     RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 16     RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#> 18     RELATED Ailuropoda Ursidae Carnivora Mammalia Chordata  VULNERABLE
#>    sp_nameMatch
#> 1         INPUT
#> 2         EXACT
#> 3         EXACT
#> 4         EXACT
#> 5         EXACT
#> 6         EXACT
#> 8         EXACT
#> 9         EXACT
#> 10        EXACT
#> 14        EXACT
#> 15        EXACT
#> 16        EXACT
#> 18        EXACT

# --- 4. Fuzzy matching for uncertain names ---
# sp_nameMatch = "VARIANT" for close but non-identical fuzzy matches
get_status("Ailuropoda melanoleuca", search = FALSE)
#>                  canonicalName    rank gbif_key
#> 2433399 Ailuropoda melanoleuca SPECIES  2433399
#> 9387176 Aeluropus melanoleucus SPECIES  9387176
#> 9379787 Ailuropus melanoleucus SPECIES  9379787
#> 7888034     Ursus melanoleucus SPECIES  7888034
#>                               scientificName gbif_status      Genus  Family
#> 2433399 Ailuropoda melanoleuca (David, 1869)    ACCEPTED Ailuropoda Ursidae
#> 9387176 Aeluropus melanoleucus (David, 1869)     SYNONYM Ailuropoda Ursidae
#> 9379787 Ailuropus melanoleucus (David, 1869)     SYNONYM Ailuropoda Ursidae
#> 7888034       Ursus melanoleucus David, 1869     SYNONYM Ailuropoda Ursidae
#>             Order    Class   Phylum IUCN_status sp_nameMatch
#> 2433399 Carnivora Mammalia Chordata  VULNERABLE        INPUT
#> 9387176 Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 9379787 Carnivora Mammalia Chordata  VULNERABLE        EXACT
#> 7888034 Carnivora Mammalia Chordata  VULNERABLE        EXACT

# --- 5. Cross-check get_status() keys against get_gbif() output ---
occ <- get_gbif("Ailuropoda melanoleuca", has_xy = TRUE, verbose = FALSE)
valid_keys    <- tax_ch$gbif_key[tax_ch$gbif_status %in% c("ACCEPTED", "CHILDREN")]
returned_keys <- unique(occ$acceptedTaxonKey)

# Keys in get_status() — backbone taxa (ACCEPTED + CHILDREN)
valid_keys
#> [1] "2433399" "9087131" "9460321"

# Keys returned by get_gbif()
returned_keys
#> [1] 9087131 2433399 9460321

# Are all occurrence keys known backbone keys?
# FALSE indicates non-backbone entries (e.g. BOLD sequences) are present
all(returned_keys %in% valid_keys)
#> [1] TRUE

# Inspect non-backbone records — typically MATERIAL_SAMPLE, excluded by default
extra_keys <- returned_keys[!returned_keys %in% valid_keys]
extra_keys
#> integer(0)
occ[occ$acceptedTaxonKey %in% extra_keys,
    c("acceptedTaxonKey", "scientificName", "basisOfRecord")]
#> [1] acceptedTaxonKey scientificName   basisOfRecord   
#> <0 rows> (or 0-length row.names)

# Note: not all valid_keys need to appear in returned_keys — some subspecies
# may have no records in GBIF at all (e.g. extinct taxa) or may have been
# removed by get_gbif() filtering options (e.g. has_xy, basis, grain)

# --- 6. Input name flagged correctly regardless of how GBIF resolves it ---
# sp_nameMatch = "INPUT" marks the row closest to the submitted name
tax_ch[tax_ch$sp_nameMatch == "INPUT", ]
#>                  canonicalName    rank gbif_key
#> 2433399 Ailuropoda melanoleuca SPECIES  2433399
#>                               scientificName gbif_status      Genus  Family
#> 2433399 Ailuropoda melanoleuca (David, 1869)    ACCEPTED Ailuropoda Ursidae
#>             Order    Class   Phylum IUCN_status sp_nameMatch
#> 2433399 Carnivora Mammalia Chordata  VULNERABLE        INPUT

# }
```
