\donttest{
# Load available ecoregions
eco.terra <- read_ecoreg(
    ecoreg_name = "eco_terra",
    save_dir = tempdir(),
    format = "sf"
)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pd <- get_gbif(sp_name = "Ailuropoda melanoleuca")

# Build a range map from occurrence points
range.panda <- get_range(
    occ_coord = obs.pd,
    ecoreg = eco.terra,
    clust_pts_outlier = 6,
    ecoreg_name = "ECO_NAME",
    format = "sf"
)
pd.test <- cv_range(
    range_object = range.panda,
    cv = "block-cv",
    nfolds = 5,
    nblocks = 2
);pd.test

}
