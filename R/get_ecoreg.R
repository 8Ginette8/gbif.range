#' Metadata for Downloadable Ecoregion Layers
#'
#' A named list describing the ecoregion layers that can be downloaded with
#' \code{get_ecoreg()}. Each element is itself a list with three character
#' fields: \code{filename} (the local save name), \code{link} (the download
#' URL), and \code{description} (a short label).
#'
#' @format A list with each element containing:
#' \describe{
#'   \item{filename}{The name to save the downloaded file as.}
#'   \item{link}{The URL to download the file.}
#'   \item{description}{A short description of the file.}
#' }
#' @return A named list of length 4, one element per available ecoregion
#' dataset, each with \code{filename}, \code{link}, and \code{description}
#' fields as described in \code{Format}.
#' @docType data
#' @seealso \code{\link{get_ecoreg}}() to download the datasets listed here,
#' and \code{\link{read_ecoreg}}() to load one directly.
#' @examples
#' ecoreg_list
#' @export
ecoreg_list <- list(
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
#' Determine the directory used to store downloaded ecoregion data.
#'
#' @param save_dir Optional character. Directory where downloaded files should be stored.
#' If \code{NULL} (the default), the directory is resolved in this order:
#' (1) the \code{gbif.range.save_dir} R option, (2) the
#' \code{GBIF_RANGE_SAVE_DIR} environment variable, (3) \code{tempdir()}.
#' No files are written outside \code{tempdir()} unless the user has
#' explicitly configured a persistent location via one of these two settings.
#' @return A character string giving the target directory.
#' @noRd
get_save_dir <- function(save_dir = NULL) {
  if (is.null(save_dir)) {
    save_dir <- getOption("gbif.range.save_dir")   # (1)
  }
  if (is.null(save_dir)) {
    env_dir <- Sys.getenv("GBIF_RANGE_SAVE_DIR", unset = "")
    if (nzchar(env_dir)) save_dir <- env_dir       # (2)
  }
  if (is.null(save_dir)) {
    save_dir <- tempdir()                          # (3)
  }
  return(save_dir)
}


#' Download Ecoregion Layers
#'
#' Download one or more ecoregion datasets listed in \code{ecoreg_list()} and
#' unpack them into a local directory.
#'
#' @param ecoreg_name Character. Use \code{"all"} to download every
#' dataset, or supply a single file name listed in \code{ecoreg_list}.
#' @param save_dir Character. Directory where the downloaded zip files and extracted
#' shapefiles should be stored. Defaults to \code{tempdir()} unless the
#' \code{gbif.range.save_dir} option or \code{GBIF_RANGE_SAVE_DIR} environment
#' variable is set, in which case that location is used instead.
#' @return \code{NULL}. Files are downloaded and unpacked for side effects.
#' @seealso \code{\link{ecoreg_list}} for the list of available datasets, and
#' \code{\link{read_ecoreg}}() to download and load a layer in one step.
#' @examples
#' \donttest{
#' # Download every ecoregion dataset listed in ecoreg_list
#' get_ecoreg(save_dir = tempdir())
#' }
#' @export
get_ecoreg <- function(ecoreg_name = "all", save_dir = NULL) {
  save_dir <- get_save_dir(save_dir)
  
  # Ensure the directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Check if ecoreg_name exists in ecoreg_list
  data_list <- ecoreg_list
  
  if (ecoreg_name != "all") {
    # If ecoreg_name does not exist, error
    va1 <- vapply(data_list, function(x) x$filename, character(1))
    if (!ecoreg_name %in% va1) {
      stop(
        paste(
          "Ecoregion", ecoreg_name, "not found. \n Available ecoregions are:",
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
    data_list <- data_list[va2 == ecoreg_name]
    message(
      "Preparing to download ecoregion ",
      ecoreg_name,
      " \n data file to: ",
      save_dir
    )
  } else {
    message(
      "Preparing to download all listed ecoregions: ",
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


#' Check for a Local Ecoregion Layer and Download It if Needed
#'
#' Check whether a target ecoregion directory exists and contains at least one
#' shapefile. If not, download the dataset with \code{get_ecoreg()}.
#'
#' @param ecoreg_name Character. File name of the ecoregion dataset to check. See
#' \code{ecoreg_list} for valid values.
#' @param save_dir Character. Directory where downloaded files should be stored.
#' Defaults to \code{tempdir()} unless the \code{gbif.range.save_dir} option or
#' \code{GBIF_RANGE_SAVE_DIR} environment variable is set, in which case that
#' location is used instead.
#' @return \code{NULL}. The function downloads data only when required.
#' @seealso \code{\link{get_ecoreg}}() to force a download, and
#' \code{\link{read_ecoreg}}() which calls this function internally before
#' loading the layer.
#' @examples
#' \donttest{
#' check_and_get_ecoreg("eco_terra", save_dir = tempdir())
#' }
#' @export
check_and_get_ecoreg <- function(ecoreg_name = "eco_terra", save_dir = NULL) {
  save_dir <- get_save_dir(save_dir)
  
  # Determine the directory to check
  check_dir <- file.path(save_dir, ecoreg_name)
  
  # Check if directory exists and contains at least one .shp file
  if (
    !(dir.exists(check_dir) &&
      length(list.files(check_dir, pattern = "\\.shp$", full.names = TRUE)) > 0)
  ) {
    message(
      ecoreg_name,
      " directory does not exist or contains no .shp files. 
            \n [", file.path(check_dir, ecoreg_name),
      "] will be created and data will be downloaded.
            \n Downloading data...")
    get_ecoreg(ecoreg_name, save_dir)
  }
}


#' Read an Ecoregion Layer
#' 
#' Load one of the packaged ecoregion datasets listed in \code{ecoreg_list()}.
#' If the requested data are not available locally, they are downloaded first.
#'
#' @param ecoreg_name Character. File name of the ecoregion dataset to load. See
#' \code{ecoreg_list} for available options.
#' @param save_dir Character. Directory where the downloaded files are stored.
#' Defaults to \code{tempdir()} unless the \code{gbif.range.save_dir} option or
#' \code{GBIF_RANGE_SAVE_DIR} environment variable is set, in which case that
#' location is used instead.
#' @param format \code{"SpatVector"} (default) or
#' \code{"sf"}.
#' @details Four datasets are currently available:
#'
#' (1) \code{eco_terra} for terrestrial species, based on The Nature
#' Conservancy version adapted from Olson et al. (2001).
#' 
#' (2) \code{eco_marine} for marine species, based on Spalding et al. (2007,
#' 2012).
#' 
#' (3) \code{eco_hd_marine}, a higher-resolution marine version.
#' 
#' (4) \code{eco_fresh} for freshwater species, based on Abell et al. (2008).
#' @return An ecoregion layer as a \code{SpatVector} or \code{sf} object,
#' depending on \code{format}.
#' @references
#' Olson, D. M., Dinerstein, E., Wikramanayake, E. D., Burgess, N. D.,
#' Powell, G. V. N., Underwood, E. C., D'Amico, J. A., Itoua, I.,
#' Strand, H. E., Morrison, J. C., Loucks, C. J., Allnutt, T. F.,
#' Ricketts, T. H., Kura, Y., Lamoreux, J. F., Wettengel, W. W., Hedao,
#' P., Kassem, K. R. 2001. Terrestrial ecoregions of the world: a new map
#' of life on Earth. BioScience 51(11):933-938.
#' \doi{10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2}
#' 
#' The Nature Conservancy (2009). Global Ecoregions, Major Habitat Types,
#' Biogeographical Realms and The Nature Conservancy Terrestrial Assessment
#' Units. GIS layers developed by The Nature Conservancy with multiple partners,
#' combined from Olson et al. (2001), Bailey 1995 and Wiken 1986. Cambridge
#' (UK): The Nature Conservancy.
#' 
#' Spalding, M. D., Fox, H. E., Allen, G. R., Davidson, N., Ferdana, Z. A.,
#' Finlayson, M., Halpern, B. S., Jorge, M. A., Lombana, A., Lourie, S. A.,
#' Martin, K. D., McManus, E., Molnar, J., Recchia, C. A., Robertson, J.
#' (2007). Marine Ecoregions of the World: A Bioregionalization of Coastal
#' and Shelf Areas. BioScience, 57(7), 573-583. \doi{10.1641/B570707}
#' 
#' Spalding, M. D., Agostini, V. N., Rice, J., & Grant, S. M. (2012).
#' Pelagic provinces of the world: a biogeographic classification of the
#' world's surface pelagic waters. Ocean & Coastal Management, 60, 19-30.
#' \doi{10.1016/j.ocecoaman.2011.12.016}
#' 
#' The Nature Conservancy (2012). Marine Ecoregions and Pelagic Provinces
#' of the World. GIS layers developed by The Nature Conservancy with multiple
#' partners, combined from Spalding et al. (2007) and Spalding et al. (2012).
#' Cambridge (UK): The Nature Conservancy.
#' 
#' Abell, R., Thieme, M. L., Revenga, C., Bryer, M., Kottelat, M.,
#' Bogutskaya, N., Coad, B., Mandrak, N., Contreras Balderas, S., Bussing,
#' W., Stiassny, M. L. J., Skelton, P., Allen, G. R., Unmack, P., Naseka,
#' A., Ng, R., Sindorf, N., Robertson, J., Armijo, E., Higgins, J. V.,
#' Heibel, T. J., Wikramanayake, E., Olson, D., Lopez, H. L., Reis, R. E.,
#' Lundberg, J. G., Sabaj Perez, M. H., Petry, P. (2008). Freshwater
#' Ecoregions of the World: A New Map of Biogeographic Units for Freshwater
#' Biodiversity Conservation. BioScience, 58(5), 403-414.
#' \doi{10.1641/B580507}
#' @seealso \code{\link{ecoreg_list}} for the list of available datasets,
#' \code{\link{get_ecoreg}}() to force a fresh download, and
#' \code{\link{check_and_get_ecoreg}}() for the underlying download check.
#' @examples
#' \donttest{
#' shp_eco_terra <- read_ecoreg("eco_terra", save_dir = tempdir())
#' terra::plot(shp_eco_terra)
#' }
#' @export
read_ecoreg <- function(ecoreg_name = "eco_terra",
                        save_dir = NULL,
                        format = "SpatVector") {
    # First function
    save_dir <- get_save_dir(save_dir)
  
    # check and if non existing get ecoreg
    check_and_get_ecoreg(ecoreg_name, save_dir)

    # Load the shapefile
    shp_path <- list.files(
     file.path(save_dir, ecoreg_name),
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
