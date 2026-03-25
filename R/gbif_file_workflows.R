### =========================================================================
### File-based GBIF workflows
### =========================================================================
#' Split a Downloaded GBIF Table into One File per Species
#'
#' Stream a large GBIF export from disk in chunks and write one occurrence file
#' per species or GBIF taxon key. The function is designed for multi-species
#' tables that are too large to load fully into memory.
#'
#' @param input_file Path to a tabular GBIF export already stored on disk.
#' @param outdir Directory where one per-species file will be written.
#' @param chunk_size Number of rows read at a time. Larger values are usually
#'   faster, whereas smaller values reduce peak memory use.
#' @param select_cols Character vector of columns to keep from the original
#'   file. The defaults retain the taxon key, species labels, and geographic
#'   coordinates needed by downstream range workflows.
#' @param sep_in Field separator used by the input file. GBIF downloads are
#'   usually tab-delimited.
#' @param sep_out Field separator used for the saved species files.
#' @param overwrite Logical. If \code{TRUE}, existing batch files created by
#'   this function in \code{outdir} are removed before writing new ones.
#' @param verbose Logical. Should progress messages be printed?
#' @details Output file names follow the pattern
#'   \code{occurrences_speciesKey_<key>_<species>.csv}. The extension is kept as
#'   \code{.csv} for convenience, even when the file remains tab-delimited.
#' @return A data frame summarizing the written files, with one row per species
#'   key and the columns \code{species_key}, \code{species_name},
#'   \code{n_records}, and \code{species_file}.
#' @seealso \code{\link{species_csvs_to_ranges}()} to process the written
#'   species files sequentially with \code{\link{get_range}()}.
#' @example inst/examples/split_gbif_by_species_help.R
#' @export
split_gbif_by_species <- function(
    input_file,
    outdir = file.path(tempdir(), "gbif_by_species"),
    chunk_size = 100000,
    select_cols = c(
      "speciesKey",
      "species",
      "scientificName",
      "decimalLongitude",
      "decimalLatitude"
    ),
    sep_in = "\t",
    sep_out = "\t",
    overwrite = FALSE,
    verbose = TRUE) {

  gbif_require_data_table("split_gbif_by_species")
  check_character_vector(input_file, "input_file")
  check_character_vector(outdir, "outdir")
  check_numeric(chunk_size, "chunk_size")
  check_character_vector(select_cols, "select_cols")
  check_character_vector(sep_in, "sep_in")
  check_character_vector(sep_out, "sep_out")
  check_logical(overwrite, "overwrite")
  check_logical(verbose, "verbose")

  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  # Prepare the output directory once so every chunk can append safely.
  gbif_prepare_output_dir(
    outdir = outdir,
    pattern = "^occurrences_speciesKey_.*\\.csv$",
    overwrite = overwrite,
    label = "species occurrence files"
  )

  # Read the header only once, then reuse the column indices while streaming
  # the file chunk by chunk.
  header <- data.table::fread(
    input_file,
    sep = sep_in,
    nrows = 0,
    showProgress = FALSE,
    fill = TRUE
  )
  all_cols <- names(header)

  required_cols <- c("speciesKey", "decimalLongitude", "decimalLatitude")
  missing_cols <- setdiff(required_cols, all_cols)
  if (length(missing_cols) > 0) {
    stop(
      "The input file is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  cols_to_read <- unique(intersect(select_cols, all_cols))
  if (length(cols_to_read) == 0) {
    stop("None of the requested 'select_cols' were found in the input file.")
  }
  select_idx <- match(cols_to_read, all_cols)

  # Keep one running summary entry per species key across all chunks.
  summary_map <- new.env(parent = emptyenv())
  rows_read <- 0L

  repeat {
    # Read the next slice of the file without loading previous chunks again.
    chunk <- gbif_fread_chunk(
      input_file = input_file,
      sep = sep_in,
      col_names = all_cols,
      select_idx = select_idx,
      rows_read = rows_read,
      chunk_size = chunk_size
    )

    if (nrow(chunk) == 0) {
      break
    }
    chunk_n <- nrow(chunk)

    chunk <- as.data.frame(chunk, stringsAsFactors = FALSE)
    # Only keyed occurrences with usable coordinates can contribute to a
    # species-level occurrence file or a downstream range.
    chunk <- chunk[!is.na(chunk$speciesKey), , drop = FALSE]
    if (nrow(chunk) > 0) {
      chunk <- chunk[
        !is.na(chunk$decimalLongitude) & !is.na(chunk$decimalLatitude),
        ,
        drop = FALSE
      ]
    }

    if (nrow(chunk) > 0) {
      chunk$speciesKey <- gbif_key_to_string(chunk$speciesKey)
      chunk$species_name <- gbif_choose_species_name(
        species = chunk[["species"]],
        scientific_name = chunk[["scientificName"]]
      )

      # Append rows species by species so the output directory becomes a set of
      # single-species files ready for get_range()-style processing.
      for (key in unique(chunk$speciesKey)) {
        idx <- which(chunk$speciesKey %in% key)
        species_name <- gbif_choose_name_from_chunk(chunk$species_name[idx], key)
        outfile <- file.path(
          outdir,
          gbif_split_file_name(key = key, species_name = species_name)
        )

        to_write <- chunk[idx, setdiff(names(chunk), "species_name"), drop = FALSE]
        gbif_write_delim(to_write, outfile, sep = sep_out)
        gbif_record_split_summary(summary_map, key, species_name, outfile, nrow(to_write))
      }
    }

    rows_read <- rows_read + chunk_n
    if (chunk_n < chunk_size) {
      break
    }

    if (verbose) {
      message("Processed ", rows_read, " rows from ", basename(input_file), ".")
    }
  }

  gbif_split_summary_df(summary_map)
}


#' Build Range Maps from Per-Species GBIF Files Saved on Disk
#'
#' Read one occurrence file per species, prepare the minimal coordinate table
#' needed by \code{\link{get_range}()}, and save one range output per species.
#'
#' @param species_dir Directory containing files created by
#'   \code{\link{split_gbif_by_species}()}.
#' @param ecoreg Ecoregion input accepted by \code{\link{get_range}()}. This
#'   can be a spatial object already loaded in memory, or a built-in ecoregion
#'   name such as \code{"eco_terra"} that will be resolved with
#'   \code{\link{read_ecoreg}()}.
#' @param ecoreg_name Name of the categorical field defining the ecoregion
#'   units. For example, use \code{"ECO_NAME"} with \code{"eco_terra"}.
#' @param outdir Directory where range files will be saved.
#' @param occ_outdir Optional directory where the minimal occurrence tables
#'   passed to \code{\link{get_range}()} will also be saved.
#' @param occ_save_as File format for minimal occurrence tables. Choose among
#'   \code{"none"}, \code{"tsv"}, or \code{"rds"}.
#' @param range_save_as File format for range outputs. Choose among
#'   \code{"rds"}, \code{"gpkg"}, or \code{"tif"}.
#' @param deduplicate Logical. Should identical longitude-latitude pairs be
#'   collapsed before range inference?
#' @param sep_in Field separator used by the per-species occurrence files.
#' @param overwrite Logical. If \code{TRUE}, existing outputs in
#'   \code{outdir}/\code{occ_outdir} are replaced.
#' @param verbose Logical. Should progress messages be printed?
#' @param ... Additional arguments passed to \code{\link{get_range}()}.
#' @details The function is intended to work seamlessly with
#'   \code{\link{split_gbif_by_species}()}. It reads each species file
#'   sequentially, adds the single-species \code{input_search} column expected
#'   by \code{get_range()}, and keeps only the coordinate columns needed for
#'   range construction.
#' @return A data frame summarizing the processed files, with one row per range
#'   and the columns \code{species_key}, \code{species_name}, \code{n_points},
#'   \code{occ_file}, and \code{range_file}.
#' @seealso \code{\link{split_gbif_by_species}()} and
#'   \code{\link{read_range_rds}()}.
#' @example inst/examples/species_csvs_to_ranges_help.R
#' @export
species_csvs_to_ranges <- function(
    species_dir,
    ecoreg,
    ecoreg_name = NULL,
    outdir = file.path(tempdir(), "gbif_ranges"),
    occ_outdir = NULL,
    occ_save_as = c("none", "tsv", "rds"),
    range_save_as = c("rds", "gpkg", "tif"),
    deduplicate = TRUE,
    sep_in = "\t",
    overwrite = FALSE,
    verbose = TRUE,
    ...) {

  gbif_require_data_table("species_csvs_to_ranges")
  check_character_vector(species_dir, "species_dir")
  check_character_vector(outdir, "outdir")
  check_logical(deduplicate, "deduplicate")
  check_character_vector(sep_in, "sep_in")
  check_logical(overwrite, "overwrite")
  check_logical(verbose, "verbose")

  occ_save_as <- match.arg(occ_save_as)
  range_save_as <- match.arg(range_save_as)

  if (!dir.exists(species_dir)) {
    stop("Species directory does not exist: ", species_dir)
  }

  species_files <- list.files(
    species_dir,
    pattern = "\\.(csv|tsv)$",
    full.names = TRUE
  )
  if (length(species_files) == 0) {
    stop("No species files found in: ", species_dir)
  }

  gbif_prepare_output_dir(
    outdir = outdir,
    pattern = "^range_speciesKey_.*\\.(rds|gpkg|tif)$",
    overwrite = overwrite,
    label = "range files"
  )

  if (!is.null(occ_outdir) && !identical(occ_save_as, "none")) {
    gbif_prepare_output_dir(
      outdir = occ_outdir,
      pattern = "^occ_min_speciesKey_.*\\.(tsv|rds)$",
      overwrite = overwrite,
      label = "minimal occurrence files"
    )
  }

  # Resolve built-in ecoregion shortcuts once, outside the per-species loop.
  ecoreg_object <- gbif_resolve_ecoreg_input(ecoreg)
  summary_list <- vector("list", length(species_files))

  for (i in seq_along(species_files)) {
    file_i <- species_files[i]
    file_meta <- gbif_parse_split_filename(file_i)

    # Read only the label and coordinate columns needed to derive a minimal
    # get_range() input for this species.
    header_i <- data.table::fread(
      file_i,
      sep = sep_in,
      nrows = 0,
      showProgress = FALSE,
      fill = TRUE
    )
    cols_to_read <- intersect(
      c(
        "speciesKey",
        "species",
        "scientificName",
        "decimalLongitude",
        "decimalLatitude"
      ),
      names(header_i)
    )

    occ_i <- data.table::fread(
      file_i,
      sep = sep_in,
      select = cols_to_read,
      showProgress = FALSE,
      fill = TRUE
    )
    occ_i <- as.data.frame(occ_i, stringsAsFactors = FALSE)

    species_name <- gbif_choose_name_from_chunk(
      candidates = c(file_meta$species_name, occ_i$species, occ_i$scientificName),
      fallback = file_meta$species_key
    )

    # Collapse the per-species file to the minimal structure expected by
    # get_range(): one focal species plus decimal longitude and latitude.
    occ_min <- gbif_prepare_occ_min(
      occ = occ_i,
      species_name = species_name,
      deduplicate = deduplicate
    )

    if (nrow(occ_min) == 0) {
      if (verbose) {
        message("Skipping ", basename(file_i), " because no valid coordinates remain.")
      }
      next
    }

    occ_file <- NA_character_
    if (!is.null(occ_outdir) && !identical(occ_save_as, "none")) {
      occ_file <- gbif_save_occ_min(
        occ_min = occ_min,
        outdir = occ_outdir,
        species_key = file_meta$species_key,
        species_name = species_name,
        save_as = occ_save_as
      )
    }

    # Delegate the actual range inference to get_range() so the batch workflow
    # stays a thin file-handling layer over the core mapping algorithm.
    range_obj <- get_range(
      occ_coord = occ_min,
      ecoreg = ecoreg_object,
      ecoreg_name = ecoreg_name,
      verbose = verbose,
      ...
    )

    range_file <- gbif_save_range_output(
      range_obj = range_obj,
      outdir = outdir,
      species_key = file_meta$species_key,
      species_name = species_name,
      save_as = range_save_as,
      overwrite = overwrite
    )

    summary_list[[i]] <- data.frame(
      species_key = file_meta$species_key,
      species_name = species_name,
      n_points = nrow(occ_min),
      occ_file = occ_file,
      range_file = range_file,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message("Saved range for ", species_name, " to ", basename(range_file), ".")
    }
  }

  gbif_bind_batch_summary(summary_list)
}


#' Read a Range File Saved by \code{species_csvs_to_ranges()}
#'
#' Convenience wrapper around \code{readRDS()} for range files saved with
#' \code{range_save_as = "rds"}.
#'
#' @param file Path to an \code{.rds} range file created by
#'   \code{\link{species_csvs_to_ranges}()}.
#' @return A list with \code{init.args} and \code{rangeOutput}, matching the
#'   structure written by the batch workflow.
#' @example inst/examples/read_range_rds_help.R
#' @export
read_range_rds <- function(file) {
  check_character_vector(file, "file")
  if (!file.exists(file)) {
    stop("Range file does not exist: ", file)
  }

  obj <- readRDS(file)
  if (!is.list(obj) || !"rangeOutput" %in% names(obj)) {
    stop("The supplied RDS file is not a gbif.range batch range file.")
  }

  if (inherits(obj$rangeOutput, "gbifPackedSpatVector")) {
    obj$rangeOutput <- terra::unwrap(obj$rangeOutput$packed)
  } else if (inherits(obj$rangeOutput, "gbifPackedSpatRaster")) {
    obj$rangeOutput <- terra::unwrap(obj$rangeOutput$packed)
  } else if (inherits(obj$rangeOutput, "PackedSpatVector")) {
    obj$rangeOutput <- terra::unwrap(obj$rangeOutput)
  } else if (inherits(obj$rangeOutput, "PackedSpatRaster")) {
    obj$rangeOutput <- terra::unwrap(obj$rangeOutput)
  }

  obj
}


#' Plot a Packed Vector Range Saved by \code{species_csvs_to_ranges()}
#'
#' @method plot gbifPackedSpatVector
#' @param x Object of class \code{gbifPackedSpatVector}.
#' @param ... Additional arguments passed to \code{terra::plot()}.
#' @keywords internal
#' @export
plot.gbifPackedSpatVector <- function(x, ...) {
  terra::plot(terra::unwrap(x$packed), ...)
}


#' Plot a Packed Raster Range Saved by \code{species_csvs_to_ranges()}
#'
#' @method plot gbifPackedSpatRaster
#' @param x Object of class \code{gbifPackedSpatRaster}.
#' @param ... Additional arguments passed to \code{terra::plot()}.
#' @keywords internal
#' @export
plot.gbifPackedSpatRaster <- function(x, ...) {
  terra::plot(terra::unwrap(x$packed), ...)
}


#' Check for an Optional \code{data.table} Dependency
#'
#' @param caller Character string naming the calling function.
#' @keywords internal
#' @noRd
gbif_require_data_table <- function(caller) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package 'data.table' is required for ", caller, "(). ",
      "Please install it with install.packages('data.table')."
    )
  }
}


#' Prepare an Output Directory for a Fresh Batch Run
#'
#' @param outdir Output directory.
#' @param pattern File pattern to clear when \code{overwrite = TRUE}.
#' @param overwrite Logical. Should existing matching files be removed?
#' @param label Human-readable output label for error messages.
#' @keywords internal
#' @noRd
gbif_prepare_output_dir <- function(outdir, pattern, overwrite, label) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  gbif_remove_existing_files(
    outdir = outdir,
    pattern = pattern,
    overwrite = overwrite,
    label = label
  )
}


#' Read One Chunk from a Delimited File with \code{data.table}
#'
#' @param input_file Input file path.
#' @param sep Field separator.
#' @param col_names Full set of input column names.
#' @param select_idx Integer indices of the columns to keep.
#' @param rows_read Number of data rows already processed.
#' @param chunk_size Number of rows to read.
#' @keywords internal
#' @noRd
gbif_fread_chunk <- function(
    input_file,
    sep,
    col_names,
    select_idx,
    rows_read,
    chunk_size) {

  empty_chunk <- function() {
    out <- as.data.frame(
      matrix(nrow = 0, ncol = length(select_idx)),
      stringsAsFactors = FALSE
    )
    names(out) <- col_names[select_idx]
    out
  }

  if (rows_read == 0L) {
    return(data.table::fread(
      input_file,
      sep = sep,
      select = select_idx,
      nrows = chunk_size,
      showProgress = FALSE,
      fill = TRUE
    ))
  }

  tryCatch(
    data.table::fread(
      input_file,
      sep = sep,
      skip = rows_read + 1L,
      header = FALSE,
      col.names = col_names[select_idx],
      select = select_idx,
      nrows = chunk_size,
      showProgress = FALSE,
      fill = TRUE
    ),
    error = function(e) {
      if (grepl("skip=", conditionMessage(e), fixed = TRUE)) {
        return(empty_chunk())
      }
      stop(e)
    }
  )
}


#' Convert GBIF Taxon Keys to Stable Character Strings
#'
#' @param x Taxon key values.
#' @keywords internal
#' @noRd
gbif_key_to_string <- function(x) {
  out <- trimws(format(x, scientific = FALSE, trim = TRUE))
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out
}


#' Sanitize Species Labels for File Names
#'
#' @param x Character vector of species labels.
#' @keywords internal
#' @noRd
gbif_safe_file_name <- function(x) {
  x <- trimws(x)
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^[:alnum:]_.-]", "", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[!nzchar(x)] <- "unknown_species"
  x
}


#' Choose a Species Label from GBIF Name Columns
#'
#' @param species Vector from the GBIF \code{species} column.
#' @param scientific_name Vector from the GBIF \code{scientificName} column.
#' @keywords internal
#' @noRd
gbif_choose_species_name <- function(species, scientific_name) {
  species <- if (is.null(species)) rep(NA_character_, length(scientific_name)) else species
  scientific_name <- if (is.null(scientific_name)) rep(NA_character_, length(species)) else scientific_name

  out <- ifelse(
    !is.na(species) & nzchar(trimws(species)),
    trimws(species),
    trimws(scientific_name)
  )
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out
}


#' Choose the Best Available Species Label from a Set of Candidates
#'
#' @param candidates Character vector of possible species labels.
#' @param fallback Fallback label if all candidates are empty.
#' @keywords internal
#' @noRd
gbif_choose_name_from_chunk <- function(candidates, fallback) {
  candidates <- unique(trimws(as.character(candidates)))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 0) {
    return(as.character(fallback))
  }
  candidates[which.max(nchar(candidates))]
}


#' Build the Standard Per-Species File Name
#'
#' @param key GBIF taxon key.
#' @param species_name Species label used in the file name.
#' @keywords internal
#' @noRd
gbif_split_file_name <- function(key, species_name) {
  paste0(
    "occurrences_speciesKey_",
    gbif_key_to_string(key),
    "_",
    gbif_safe_file_name(species_name),
    ".csv"
  )
}


#' Parse Metadata Back from a Per-Species File Name
#'
#' @param path File path created by \code{split_gbif_by_species()}.
#' @keywords internal
#' @noRd
gbif_parse_split_filename <- function(path) {
  file_stub <- tools::file_path_sans_ext(basename(path))
  hit <- regexec("^occurrences_speciesKey_([^_]+)_(.*)$", file_stub)
  parts <- regmatches(file_stub, hit)[[1]]

  if (length(parts) >= 3) {
    list(
      species_key = parts[2],
      species_name = gsub("_", " ", parts[3], fixed = TRUE)
    )
  } else {
    list(
      species_key = NA_character_,
      species_name = NA_character_
    )
  }
}


#' Append Delimited Output Without Repeating the Header
#'
#' @param x Data frame to write.
#' @param file Output file path.
#' @param sep Field separator.
#' @keywords internal
#' @noRd
gbif_write_delim <- function(x, file, sep) {
  utils::write.table(
    x,
    file = file,
    sep = sep,
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file),
    quote = FALSE,
    na = ""
  )
}


#' Remove Existing Batch Outputs if Requested
#'
#' @param outdir Output directory.
#' @param pattern File-name pattern to match.
#' @param overwrite Logical. Should existing files be removed?
#' @param label Human-readable file label used in error messages.
#' @keywords internal
#' @noRd
gbif_remove_existing_files <- function(outdir, pattern, overwrite, label) {
  existing <- list.files(outdir, pattern = pattern, full.names = TRUE)
  if (length(existing) == 0) {
    return(invisible(NULL))
  }
  if (!overwrite) {
    stop(
      "Existing ", label, " were found in ", outdir,
      ". Set overwrite = TRUE to replace them."
    )
  }
  file.remove(existing)
  invisible(NULL)
}


#' Record Split-File Metadata in an Environment
#'
#' @param summary_map Environment used to accumulate metadata.
#' @param key Species key.
#' @param species_name Species label.
#' @param outfile Output file path.
#' @param n_records Number of rows appended in the current chunk.
#' @keywords internal
#' @noRd
gbif_record_split_summary <- function(summary_map, key, species_name, outfile, n_records) {
  if (!exists(key, envir = summary_map, inherits = FALSE)) {
    assign(
      key,
      list(
        species_key = key,
        species_name = species_name,
        n_records = n_records,
        species_file = outfile
      ),
      envir = summary_map
    )
  } else {
    current <- get(key, envir = summary_map, inherits = FALSE)
    current$n_records <- current$n_records + n_records
    if (nchar(species_name) > nchar(current$species_name)) {
      current$species_name <- species_name
    }
    assign(key, current, envir = summary_map)
  }
}


#' Convert an Accumulator Environment to a Summary Data Frame
#'
#' @param summary_map Environment used to store split-file metadata.
#' @keywords internal
#' @noRd
gbif_split_summary_df <- function(summary_map) {
  keys <- ls(summary_map, all.names = TRUE)
  if (length(keys) == 0) {
    return(data.frame(
      species_key = character(),
      species_name = character(),
      n_records = integer(),
      species_file = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- lapply(keys, function(key) {
    as.data.frame(get(key, envir = summary_map, inherits = FALSE), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$species_name), , drop = FALSE]
}


#' Resolve the Ecoregion Input for Batch Range Building
#'
#' @param ecoreg Ecoregion object or single built-in name.
#' @keywords internal
#' @noRd
gbif_resolve_ecoreg_input <- function(ecoreg) {
  if (inherits(ecoreg, c("SpatVector", "sf", "SpatialPolygonsDataFrame", "SpatialPolygons"))) {
    return(ecoreg)
  }

  if (is.character(ecoreg) && length(ecoreg) == 1) {
    if (file.exists(ecoreg)) {
      return(terra::vect(ecoreg))
    }
    return(read_ecoreg(ecoreg_name = ecoreg, format = "SpatVector"))
  }

  stop(
    "'ecoreg' must be a spatial object accepted by get_range(), ",
    "a built-in ecoregion name such as 'eco_terra', or a readable file path."
  )
}


#' Build the Minimal Occurrence Table Needed by \code{get_range()}
#'
#' @param occ Data frame of per-species occurrences.
#' @param species_name Single species label assigned to every row.
#' @param deduplicate Logical. Should identical coordinates be collapsed?
#' @keywords internal
#' @noRd
gbif_prepare_occ_min <- function(occ, species_name, deduplicate) {
  keep <- c("decimalLongitude", "decimalLatitude")
  occ <- occ[stats::complete.cases(occ[, keep, drop = FALSE]), keep, drop = FALSE]
  if (deduplicate && nrow(occ) > 0) {
    occ <- unique(occ)
  }
  occ$input_search <- species_name
  occ <- occ[, c("input_search", "decimalLongitude", "decimalLatitude"), drop = FALSE]
  rownames(occ) <- NULL
  occ
}


#' Bind Non-Empty Per-Species Summaries into One Data Frame
#'
#' @param summary_list List of optional summary rows.
#' @keywords internal
#' @noRd
gbif_bind_batch_summary <- function(summary_list) {
  keep <- !vapply(summary_list, is.null, logical(1))
  if (!any(keep)) {
    return(data.frame(
      species_key = character(),
      species_name = character(),
      n_points = integer(),
      occ_file = character(),
      range_file = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, summary_list[keep])
}


#' Save the Minimal Occurrence Table Used for Range Inference
#'
#' @param occ_min Minimal occurrence data frame.
#' @param outdir Output directory.
#' @param species_key Species key.
#' @param species_name Species label.
#' @param save_as Output format.
#' @keywords internal
#' @noRd
gbif_save_occ_min <- function(occ_min, outdir, species_key, species_name, save_as) {
  stub <- paste0(
    "occ_min_speciesKey_",
    gbif_key_to_string(species_key),
    "_",
    gbif_safe_file_name(species_name)
  )

  if (identical(save_as, "tsv")) {
    outfile <- file.path(outdir, paste0(stub, ".tsv"))
    utils::write.table(
      occ_min,
      file = outfile,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  } else {
    outfile <- file.path(outdir, paste0(stub, ".rds"))
    saveRDS(occ_min, outfile)
  }

  outfile
}


#' Make a Range Output Safe for \code{saveRDS()}
#'
#' @param x Range output object.
#' @keywords internal
#' @noRd
gbif_make_rds_safe <- function(x) {
  if (inherits(x, c("SpatVector", "SpatRaster"))) {
    # terra::wrap() stores a portable representation that can be written to RDS.
    packed <- terra::wrap(x)
    if (inherits(x, "SpatVector")) {
      return(structure(list(packed = packed), class = "gbifPackedSpatVector"))
    }
    return(structure(list(packed = packed), class = "gbifPackedSpatRaster"))
  }
  x
}


#' Save One Range Output to Disk
#'
#' @param range_obj Result returned by \code{get_range()}.
#' @param outdir Output directory.
#' @param species_key Species key.
#' @param species_name Species label.
#' @param save_as Output format.
#' @param overwrite Logical. Passed to the spatial writer when relevant.
#' @keywords internal
#' @noRd
gbif_save_range_output <- function(
    range_obj,
    outdir,
    species_key,
    species_name,
    save_as,
    overwrite) {

  stub <- paste0(
    "range_speciesKey_",
    gbif_key_to_string(species_key),
    "_",
    gbif_safe_file_name(species_name)
  )

  if (identical(save_as, "rds")) {
    outfile <- file.path(outdir, paste0(stub, ".rds"))
    # Store a simple list rather than the live reference object so the batch
    # files can be read back safely in a new R session.
    safe_obj <- list(
      init.args = range_obj$init.args,
      rangeOutput = gbif_make_rds_safe(range_obj$rangeOutput)
    )
    saveRDS(safe_obj, outfile)
    return(outfile)
  }

  if (identical(save_as, "gpkg")) {
    outfile <- file.path(outdir, paste0(stub, ".gpkg"))
    if (inherits(range_obj$rangeOutput, "sf")) {
      sf::st_write(range_obj$rangeOutput, outfile, delete_dsn = overwrite, quiet = TRUE)
    } else if (inherits(range_obj$rangeOutput, "SpatVector")) {
      terra::writeVector(range_obj$rangeOutput, outfile, overwrite = overwrite)
    } else {
      stop("range_save_as = 'gpkg' requires a vector output from get_range().")
    }
    return(outfile)
  }

  outfile <- file.path(outdir, paste0(stub, ".tif"))
  if (!inherits(range_obj$rangeOutput, "SpatRaster")) {
    stop("range_save_as = 'tif' requires get_range(raster = TRUE).")
  }
  terra::writeRaster(range_obj$rangeOutput, outfile, overwrite = overwrite)
  outfile
}
