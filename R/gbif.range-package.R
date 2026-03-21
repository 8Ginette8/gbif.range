#' gbif.range: Tools for GBIF Retrieval and Ecoregion-Based Species Range Mapping
#'
#' \code{gbif.range} provides a workflow to retrieve occurrence records from
#' GBIF, clean and filter them for spatial analyses, assemble or download
#' ecoregion layers, generate ecologically informed species range maps, and
#' evaluate the resulting maps against independent data.
#'
#' The package is designed around a typical workflow: (1) inspect or download
#' GBIF records with \code{get_gbif_count()} and \code{get_gbif()},
#' (2) retrieve taxonomic information with \code{get_status()},
#' (3) load packaged ecoregions with \code{read_ecoreg()} or create custom ones
#' with \code{make_ecoreg()}, (4) build range maps with \code{get_range()},
#' and (5) evaluate those maps with \code{evaluate_range()} or
#' \code{cv_range()}.
#'
#' Additional helpers include \code{get_doi()} for creating GBIF-derived
#' dataset DOIs, \code{obs_filter()} for grid-based thinning, and
#' \code{make_tiles()} for splitting study extents into GBIF-ready polygons.
#'
#' @name gbif.range
#' @aliases gbif.range gbif.range-package package-gbif.range
#' @docType package
#' @example inst/examples/gbif.range-package_help.R
#' @keywords package
"_PACKAGE"
