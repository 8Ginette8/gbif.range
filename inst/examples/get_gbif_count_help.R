\dontrun{
# Get number of observations with default filters
obs.pt <- get_gbif_count(
	sp_name = "Ailuropoda melanoleuca",
	has_xy = TRUE,
	spatial_issue = FALSE,
	geo = NULL
)

# Get the total number of observations
obs.pt <- get_gbif_count(
	sp_name = "Ailuropoda melanoleuca",
	has_xy = NULL,
	spatial_issue = NULL,
	geo = NULL
)

# Example of setting 'geo' (still considered as global)
obs.pt <- get_gbif_count(
	sp_name = "Ailuropoda melanoleuca",
	has_xy = NULL,
	spatial_issue = NULL,
	geo = terra::ext()
)

# Example of fuzzy matching when search is set to FALSE
obs.pt <- get_gbif_count(
	sp_name = "Ailuropoda melanolca",
	search = FALSE
)

# Example on the European Alps
shp.lonlat <- terra::vect(
    paste0(
        system.file(package = "gbif.range"),
        "/extdata/shp_lonlat.shp"
    )
)
obs.pt <- get_gbif_count(
	sp_name = "Arctostaphylos alpinus",
	has_xy = TRUE,
	spatial_issue = FALSE,
	geo = shp.lonlat
)

}
