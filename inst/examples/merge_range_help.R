\donttest{
# Load available ecoregions
eco_terra <- read_ecoreg(ecoreg_name = "eco_terra", save_dir = tempdir())

# Worldwide observations of the giant panda from GBIF
obs_am <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Guard against an unavailable remote source
if (gbif_have(eco_terra, obs_am)) {

# Range map: one polygon per occupied ecoregion
rng <- get_range(occ_coord = obs_am, ecoreg = eco_terra,
                 ecoreg_name = "ECO_NAME")
nrow(rng$rangeOutput)
sum(terra::expanse(rng$rangeOutput, unit = "km"))

# Dissolve into a single polygon. The area is unchanged, which shows the
# ecoregion polygons did not overlap.
merged <- merge_range(rng)
nrow(merged)
sum(terra::expanse(merged, unit = "km"))

}
}
