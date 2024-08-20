\name{evaluateRM}
\alias{evaluateRM}
\title{Evaluate the Performance of Ecoregional Range Maps (ecoRM)}
\description{
  This function evaluates the performance of Ecoregional Range Maps (ecoRM) by comparing them with validation data, which can be provided as either raster files (TIF) or shapefiles (SHP). The evaluation includes calculating sensitivity and precision metrics for each species' ecological model.
}
\usage{
evaluateRM(root.dir = NULL,
           valData.dir = NULL,
           ecoRM.dir = NULL,
           valData.type = NULL,
           verbose = TRUE,
           print.map = TRUE,
           mask = NULL,
           res.fact = NULL)
}
\arguments{
  \item{root.dir}{\code{character}. The root directory where the subdirectories for the ecoregional range map and validation data are located.}
  \item{valData.dir}{\code{character}. The subdirectory within \code{root.dir} that contains the validation data files (either shapefiles or raster files).}
  \item{ecoRM.dir}{\code{character}. The subdirectory within \code{root.dir} that contains the ecological raster models (TIFF files).}
  \item{valData.type}{\code{character}. The type of the validation data. Should be either \code{"TIFF"} for raster files or \code{"SHP"} for shapefiles.}
  \item{verbose}{\code{logical}. If \code{TRUE}, the function prints progress and evaluation summaries to the console.}
  \item{print.map}{\code{logical}. If \code{TRUE}, the function generates and saves overlay maps showing the comparison between the ecoregional range map and validation data.}
  \item{mask}{\code{RasterLayer}. Optional. A raster layer used to mask or limit the area of the analysis. If \code{NULL}, the full extent of the ecological model and validation data is used.}
  \item{res.fact}{\code{numeric}. Optional. A resolution factor to resample the validation data or mask. If \code{NULL}, no resampling is performed.}
}
\value{
  A list containing:
  \item{df_eval}{A data frame summarizing the evaluation metrics for each species, including the number of true and false presences, sensitivity, and precision.}
  \item{overlay_list}{A list of rasters representing the overlay between the ecoregional range map and validation data for each species.}
}
\details{
  The function first loads the ecoregional range map (ecoRM) and validation data, ensuring they have the same coordinate reference system (CRS) and spatial extent. It then overlays the ecoregional range map with the validation data and calculates the following evaluation metrics:
  
  \itemize{
    \item \code{nbr_pres}: The total number of presence points (either in ecoRM or validation data).
    \item \code{nbr_pres_ecoRM}: The number of presence points in the ecological model.
    \item \code{nbr_true_pres}: The number of true presence points (present in both ecoRM and validation data).
    \item \code{nbr_false_pres}: The number of false presence points (present in ecoRM but not in validation data).
    \item \code{Sen_ecoRM}: Sensitivity of the ecological model (true presences / total presences).
    \item \code{Prec_ecoRM}: Precision of the ecological model (true presences / (true presences + false presences)).
  }
  For further details on the evaluation approach and performance test of ecoregional ranges against other types of range data including species distribution models and expert range maps can be found in Pinkert et al. (2023). "The potential of ecoregional range maps for boosting taxonomic coverage in ecology and conservation"  Ecography, \emph{12}, \bold{e06794}. \href{https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/ecog.06794}
  
  If \code{print.map} is set to \code{TRUE}, the function generates PDF maps showing the overlay between the ecoregional range map and validation data for each species, and saves them in an "Output" folder within the root directory.
}
\examples{
\dontrun{
out <- evaluateRM(
  root.dir = "~/2._Andere_Projekte/EcoRM-package/Analyses",
  valData.dir = "Test_data/SHP_val_data",
  ecoRM.dir = "Test_data/EcoRM_data",
  valData.type = "SHP", # either TIF or SHP
  verbose = TRUE,
  print.map = TRUE,
  mask = NULL, # raster layer to limit analysis area
  res.fact = NULL # resolution factor for resampling
)
}
}
\seealso{
  \code{\link[sf]{st_read}}, \code{\link[raster]{raster}}, \code{\link[dplyr]{summarize}}
}
\author{
  Stefan Pinkert
}
