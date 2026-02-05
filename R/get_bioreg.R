#' Bioreg data Files
#'
#' A list of data files including download links, filenames, and descriptions.
#'
#' @format A list with each element containing:
#' \describe{
#'   \item{filename}{The name to save the downloaded file as.}
#'   \item{link}{The URL to download the file.}
#'   \item{description}{A short description of the file.}
#' }
#' @examples
#' \dontrun{
#' bioreg_list
#' }
#' @export
bioreg_list <- list(
  list(
    filename = "eco_terra",
    link = paste0(
      "https://www.arcgis.com/sharing/rest/content/items/",
      "b1636d640ede4d6ca8f5e369f2dc368b/data"
      ),
    description = "Terrestrial Ecoregions of the World"
  ),
  list(
    filename = "eco_fresh",
    link = paste0(
      "https://www.arcgis.com/sharing/rest/content/items/",
      "41aa8662254f43699e792fa194d7b3bb/data"
      ),
    description = "Freshwater Ecoregions of the World"
  ),
  list(
    filename = "eco_marine",
    link = paste0(
      "https://www.arcgis.com/sharing/rest/content/items/",
      "903c3ae05b264c00a3b5e58a4561b7e6/data"
      ),
    description = "Marine Ecoregions of the World"
  ),
  list(
    filename = "eco_hd_marine",
    link = paste0(
      "https://datadownload-production.s3.amazonaws.com/",
      "WCMC036_MEOW_PPOW_2007_2012_v1.zip"
      ),
    description = paste(
      "High Definition Marine Ecoregions and Pelagic",
      "Provinces of the World (2007, 2012)"
      )
  )
)


#' Get Save Directory
#'
#' Determines the save directory within the package structure.
#'
#' @param save_dir Optional. The directory to save the downloaded files.
#' Defaults to "inst/extdata/downloads" within the package structure.
#' @return A string representing the save directory.
#' @noRd
get_save_dir <- function(save_dir = NULL) {
  if (is.null(save_dir)) {
    save_dir <- system.file("extdata/downloads", package = "gbif.range")
    if (save_dir == "") {
      save_dir <- file.path(getwd(), "inst", "extdata", "downloads")
    }
  }
  return(save_dir)
}


#' Download Bioreg Data Files
#'
#' Downloads data files from a provided list of links and saves them
#' to a specified directory.
#'
#' @param bioreg_name Character. "all" or the filename of the desired bioregion
#' to be downloaded. See \code{bioreg_list} for available options. This
#' contains a list of data file information with elements containing `link`,
#' `filename`, and `description`.
#' @param save_dir The directory to save the downloaded files. Defaults to
#' "inst/extdata/downloads" within the package structure.
#' @return \code{NULL}. Downloads the files to the specified directory.
#' @examples
#' \dontrun{
#' # download all bioregions available in bioreg_list
#' get_bioreg()
#' }
#' @export
get_bioreg <- function(bioreg_name = "all", save_dir = NULL) {
  save_dir <- get_save_dir(save_dir)
  
  # Ensure the directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Check if bioreg_name exists in bioreg_list
  data_list <- bioreg_list
  
  if (bioreg_name != "all") {
    # If bioreg_name does not exist, error
    va1 <- vapply(data_list, function(x) x$filename, character(1))
    if (!bioreg_name %in% va1) {
      stop(
        paste(
          "Bioregion", bioreg_name, "not found. \n Available bioregions are:",
          paste(
            vapply(
              data_list,
              function(x) x$filename, character(1)
            )
          ),
        collapse = " ; "
        )
      )
    }
    va2 <- vapply(data_list, function(x) x$filename, character(1))
    data_list <- data_list[va2 == bioreg_name]
    message(
      "Preparing to download bioregion ",
      bioreg_name,
      " \n data file to: ",
      save_dir
    )
  } else {
    message(
      "Preparing to download all listed bioregions: ",
      paste(
        vapply(
          data_list,
          function(x) x$filename,
          character(1)
        ),
        collapse = " ; "
      ),
      " \n data files to: ",
      save_dir
    )
  }
  
  # Download each file
  for (data in data_list) {
    # data <- data_list[[1]]
    url <- data$link
    destfile <- file.path(save_dir, paste0(data$filename, ".zip"))
    
    # Download the file
    tryCatch({
      utils::download.file(url, destfile, mode = "wb")
      message(paste("Downloaded:", data$filename))
      if (!is.null(data$description)) {
        message(paste("Description:", data$description))
      }
    }, error = function(e) {
      message(
        paste("Failed to download:", data$filename, "\nError:", e$message)
      )
    })
    
    # Unzip to folder named after filename
    zip::unzip(destfile, exdir = file.path(save_dir, data$filename))
    # Remove the zip file
    file.remove(destfile)
    message(
      paste(
        "Unzipped:", data$filename, "\n saved to:",
        file.path(save_dir, data$filename), "\n removed: ",
        paste0(data$filename, ".zip")
      )
    )
  }
}


#' Check and Download Bioreg Data Files
#'
#' Checks if a directory exists and contains at least one .shp file.
#' If not, calls get_bioreg to download the data.
#'
#' @param bioreg_name Filename of the desired bioregion to be checked and
#' downloaded if necessary. See \code{bioreg_list} for available options.
#' @param save_dir The directory to save the downloaded files. Defaults to
#' "inst/extdata/downloads" within the package structure.
#' @return \code{NULL}. Checks and downloads the files to the specified
#' directory if necessary.
#' @examples
#' \dontrun{
#' check_and_get_bioreg("eco_terra")
#' }
#' @export
check_and_get_bioreg <- function(bioreg_name = "eco_terra", save_dir = NULL) {
  save_dir <- get_save_dir(save_dir)
  
  # Determine the directory to check
  check_dir <- file.path(save_dir, bioreg_name)
  
  # Check if directory exists and contains at least one .shp file
  if (
    !(dir.exists(check_dir) &&
      length(list.files(check_dir, pattern = "\\.shp$", full.names = TRUE)) > 0)
  ) {
    message(
      bioreg_name,
      " directory does not exist or contains no .shp files. 
            \n [", file.path(check_dir, bioreg_name),
      "] will be created and data will be downloaded.
            \n Downloading data...")
    get_bioreg(bioreg_name, save_dir)
  }
}


#' Load Bioreg Data Files
#' 
#' Loads a shapefile based on the provided bioregion name and the
#' \code{bioreg_list}.
#' @param bioreg_name Filename of the desired bioregion to be loaded.
#' See \code{bioreg_list} for available options.
#' @param save_dir The directory to save the downloaded files. Defaults to
#' "inst/extdata/downloads" within the package structure.
#' @param format Character. "\code{sf}" or "\code{SpatVector}" class for
#' layer output. Default is the \code{SpatVector} class from the \code{terra}
#' package.
#' @details Four shapefiles can be downloaded with the library:
#'
#' (1) 'eco_terra' (for terrestrial species, Nature conservancy version adapted
#' from Olson & al. 2001)
#' 
#' (2) 'eco_marine' (for marine species, Spalding & al. 2007, 2012).
#' 
#' (3) 'eco_hd_marine' (higher resolution).
#' 
#' (4) 'eco_fresh' (for freshwater species, Abell & al. 2008).
#' @return \code{SpatVector} object representing the ecoregion shapefile.
#' @references
#' Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
#' Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I.,
#' Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F.,
#' Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao,
#' P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map
#' of life on Earth. BioScience 51(11):933-938. doi: 10.1641/0006-3568(2001)051
#' 
#' The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
#' Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
#' Units. GIS layers developed by The Nature Conservancy with multiple partners,
#' combined from Olson et al. (2001), Bailey 1995 and Wiken 1986. Cambridge
#' (UK): The Nature Conservancy.
#' 
#' Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A.
#' Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana,
#' Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer Molnar, Cheri
#' A. Recchia, James Robertson, Marine Ecoregions of the World: A
#' Bioregionalization of Coastal and Shelf Areas, BioScience, Volume 57,
#' Issue 7, July 2007, Pages 573–583. doi: 10.1641/B570707
#' 
#' Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
#' Pelagic provinces of the world: a biogeographic classification of the
#' world’s surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
#' doi: 10.1016/j.ocecoaman.2011.12.016
#' 
#' The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
#' of the World. GIS layers developed by The Nature Conservancy with multiple
#' partners, combined from Spalding et al. (2007) and Spalding et al. (2012).
#' Cambridge (UK): The Nature Conservancy.
#' 
#' Robin Abell, Michele L. Thieme, Carmen Revenga, Mark Bryer, Maurice
#' Kottelat, Nina Bogutskaya, Brian Coad, Nick Mandrak, Salvador Contreras
#' Balderas, William Bussing, Melanie L. J. Stiassny, Paul Skelton, Gerald R.
#' Allen, Peter Unmack, Alexander Naseka, Rebecca Ng, Nikolai Sindorf, James
#' Robertson, Eric Armijo, Jonathan V. Higgins, Thomas J. Heibel, Eric
#' Wikramanayake, David Olson, Hugo L. López, Roberto E. Reis, John G.
#' Lundberg, Mark H. Sabaj Pérez, Paulo Petry, Freshwater Ecoregions of
#' the World: A New Map of Biogeographic Units for Freshwater Biodiversity
#' Conservation, BioScience, Volume 58, Issue 5, May 2008, Pages 403–414.
#' doi: 10.1641/B580507
#' @export
#' @examples
#' \dontrun{
#' shp_eco_terra <- read_bioreg("eco_terra")
#' plot(shp_eco_terra)
#' }
read_bioreg <- function(bioreg_name = "eco_terra",
                        save_dir = NULL,
                        format = "SpatVector"){
    # First function
    save_dir <- get_save_dir(save_dir)
  
    # check and if non existing get bioreg
    check_and_get_bioreg(bioreg_name, save_dir)

    # Load the shapefile
    shp_path <- list.files(
     file.path(save_dir, bioreg_name),
     pattern = "\\.shp$",
     full.names = TRUE
    )
    shp <- terra::vect(shp_path)

    # SHP case
    if (format == "sf") {
     shp <- sf::st_as_sf(shp)
    }
    
    return(shp)
}
