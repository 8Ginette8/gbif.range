\dontrun{
# Inspect how GBIF resolves the queried name before downloading records
# Get the taxonomy and IUCN status of Cypripedium calceolus from GBIF
tax <- get_status("Cypripedium calceolus", all = FALSE)

# Keep the accepted name and synonyms used by get_gbif()
tax[  tax$gbif_status %in% c("ACCEPTED", "SYNONYM"),
  c("scientificName", "gbif_status", "gbif_key", "sp_nameMatch", "IUCN_status")]

# The ACCEPTED row identifies the GBIF taxon concept used by get_gbif()
accepted_name <- tax$scientificName[tax$gbif_status == "ACCEPTED"][1]
accepted_key <- tax$gbif_key[tax$gbif_status == "ACCEPTED"][1]

# Download occurrences for the same input name
occ <- get_gbif(
  "Cypripedium calceolus",
  has_xy = TRUE,
  verbose = FALSE
)

# Inspect the returned record-level and accepted taxonomy side by side
head(occ[
  ,
  c(
    "input_search",
    "scientificName",
    "acceptedScientificName",
    "taxonomicStatus",
    "acceptedTaxonKey"
  )
])

# The harmonized taxon key in the occurrences should match get_status()
unique(occ$acceptedTaxonKey)
accepted_key

# Use all = TRUE only if you also want children and related names for review
get_status("Cypripedium calceolus", all = TRUE)
}
