\dontrun{
# --- 1. Default: accepted name + synonyms only ---
tax <- get_status("Panthera tigris")
tax

# --- 2. Include infra-specific taxa (subspecies, varieties) ---
# These are the same keys retrieved by get_gbif()
tax_ch <- get_status("Panthera tigris", level = "children")
tax_ch

# --- 3. Include alternative name string representations (inspection only) ---
# RELATED rows are not used by get_gbif()
tax_rel <- get_status("Panthera tigris", level = "all")
tax_rel

# --- 4. Cross-check get_status() keys against get_gbif() output ---
occ <- get_gbif("Panthera tigris", has_xy = TRUE, verbose = FALSE)
valid_keys    <- tax_ch$gbif_key[tax_ch$gbif_status %in% c("ACCEPTED", "CHILDREN")]
returned_keys <- unique(occ$acceptedTaxonKey)

# Keys in get_status() — backbone taxa (ACCEPTED + CHILDREN)
valid_keys

# Keys returned by get_gbif()
returned_keys

# Are all occurrence keys known backbone keys?
# FALSE indicates non-backbone entries (e.g. BOLD sequences) are present
all(returned_keys %in% valid_keys)

# Inspect non-backbone records — typically MATERIAL_SAMPLE, excluded by default
extra_keys <- returned_keys[!returned_keys %in% valid_keys]
extra_keys
occ[occ$acceptedTaxonKey %in% extra_keys,
    c("acceptedTaxonKey", "scientificName", "basisOfRecord")]

# Note: not all valid_keys need to appear in returned_keys — some subspecies
# may have no records in GBIF at all (e.g. extinct taxa) or may have been
# removed by get_gbif() filtering options (e.g. has_xy, basis, grain)

# --- 5. Input name flagged correctly regardless of how GBIF resolves it ---
# sp_nameMatch = "INPUT" marks the row closest to the submitted name
tax_ch[tax_ch$sp_nameMatch == "INPUT", ]

# --- 6. Fuzzy matching for uncertain names ---
# sp_nameMatch = "VARIANT" for close but non-identical fuzzy matches
get_status("Panthera tigri", search = FALSE)
}
