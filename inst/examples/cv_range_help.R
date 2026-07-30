\donttest{
# Load available ecoregions
eco_terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir(),
    format = "sf"
)

# First download the worldwide observations of the panda from GBIF
obs_am <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Both calls above depend on remote services and return NULL or an empty
# table if those are unavailable, so guard the rest of the example
if (gbif_have(eco_terra, obs_am)) {

# Build a range map from occurrence points
range_panda <- get_range(
    occ_coord = obs_am,
    ecoreg = eco_terra,
    ecoreg_name = "ECO_NAME",
    format = "sf"
)
am_test <- cv_range(
    range_object = range_panda,
    cv = "block-cv",
    nfolds = 5,
    nblocks = 2
);am_test

}
}
