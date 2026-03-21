desc <- read.dcf("DESCRIPTION")
collate_raw <- desc[1, "Collate"]

if (is.na(collate_raw) || !nzchar(collate_raw)) {
  stop("DESCRIPTION is missing a Collate field.", call. = FALSE)
}

collate <- scan(text = collate_raw, what = "character", quiet = TRUE)
collate <- gsub("^'|'$", "", collate)

r_files <- basename(list.files("R", pattern = "[.]R$", full.names = TRUE))

missing_from_collate <- setdiff(r_files, collate)
missing_from_r <- setdiff(collate, r_files)
duplicated_collate <- unique(collate[duplicated(collate)])

if (length(missing_from_collate) > 0 ||
    length(missing_from_r) > 0 ||
    length(duplicated_collate) > 0) {
  msg <- c("DESCRIPTION Collate does not match the contents of R/.")

  if (length(missing_from_collate) > 0) {
    msg <- c(
      msg,
      paste(
        "Files in R/ missing from Collate:",
        paste(sort(missing_from_collate), collapse = ", ")
      )
    )
  }

  if (length(missing_from_r) > 0) {
    msg <- c(
      msg,
      paste(
        "Files listed in Collate but missing from R/:",
        paste(sort(missing_from_r), collapse = ", ")
      )
    )
  }

  if (length(duplicated_collate) > 0) {
    msg <- c(
      msg,
      paste(
        "Duplicate Collate entries:",
        paste(sort(duplicated_collate), collapse = ", ")
      )
    )
  }

  stop(paste(msg, collapse = "\n"), call. = FALSE)
}

cat("Collate OK\n")
