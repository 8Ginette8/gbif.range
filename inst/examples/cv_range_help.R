\dontrun{
# Load available ecoregions
eco_terra <- read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

# First download the worldwide observations of Panthera tigris from GBIF
occ <- get_gbif("Panthera tigris",
          basis = c("OBSERVATION", 
                  "HUMAN_OBSERVATION", 
                  "MACHINE_OBSERVATION"))

# Make range from occurance points
range <- get_range(occ, eco_terra,"ECO_NAME")

cv.range <- cv_range(range_object = range, cv = 'block-cv', nfolds = 5, nblocks = 2)
}