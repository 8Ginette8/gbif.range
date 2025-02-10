\dontrun{
# Load available ecoregions
eco.terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# First download the worldwide observations of Panthera tigris from GBIF
obs.pt <- get_gbif(sp_name = "Panthera tigris")

# Make range from occurance points
range.tiger <- get_range(occ_coord = obs.pt,
                         bioreg = eco.terra,
                         bioreg_name = "ECO_NAME",
                         degrees_outlier = 6,
                         clustered_points_outlier = 4)

cv_range(range_object = range.tiger,
         cv = 'block-cv',
         nfolds = 5,
         nblocks = 2)
}
