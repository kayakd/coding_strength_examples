#' fetch_geo_supp_file
#'
#' @description
#' Downloads and optionally loads supplementary matrix files (e.g., .csv, .tsv, .txt) from a given NCBI GEO Series (GSE) ID
#' using the GEOquery Bioconductor package. Interactive file selection is supported.
#'
#' The function assigns the selected files as a named list to the global environment
#' under a variable named `<GSE_ID>_data`.
#'
#' @param gse_id Character. A GEO Series ID (e.g., "GSE123456"). If NULL, the user is prompted to input one interactively.
#' @param read_to_memory Logical. If TRUE (default), selected files will be read into R memory and stored as a named list.
#'
#' @return
#' If `read_to_memory = TRUE`, a list of data frames is created and assigned into the global environment.
#' If `read_to_memory = FALSE`, files are downloaded only.
#'
#' @examples
#' fetch_geo_supp_file("GSE197937")
#'
#' @author
#' Created by Koray Doğan Kaya, 2020
#' Updated by Koray Doğan Kaya, 2025

fetch_geo_supp_file <- function(gse_id = NULL, read_to_memory = TRUE) {
  # Ensure GEOquery is available
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Please install the GEOquery package first: BiocManager::install('GEOquery')")
  }
  library(GEOquery)

  # Ask user for GSE ID if not provided
  if (is.null(gse_id)) {
    gse_id <- readline(prompt = "Enter GEO Series ID (e.g., GSE123456): ")
  }

  # Set base directory for downloaded files
  base_dir <- file.path("GEO_data", gse_id)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

  # Download supplementary files
  message("Fetching available supplementary files...")
  GEOquery::getGEOSuppFiles(gse_id, baseDir = "GEO_data")

  # Get list of all downloaded files
  supp_files <- list.files(base_dir, full.names = TRUE)
  short_names <- basename(supp_files)

  if (length(supp_files) == 0) {
    stop("No supplementary files found for this GSE.")
  }

  # List available files for user selection
  cat("\nAvailable supplementary files:\n")
  for (i in seq_along(short_names)) {
    cat(sprintf("[%d] %s\n", i, short_names[i]))
  }

  # Prompt for file selection
  selection <- readline(prompt = "Enter file number(s) to load (comma-separated): ")
  indices <- as.integer(strsplit(selection, ",")[[1]])
  selected_files <- supp_files[indices]

  # Initialize list to store parsed data
  loaded_data <- list()

  for (file_path in selected_files) {
    file_name <- basename(file_path)
    cat(sprintf("\nReading: %s\n", file_name))

    if (read_to_memory) {
      if (grepl("\\.gz$", file_name)) {
        con <- gzfile(file_path, open = "rt")
        if (grepl("\\.csv", file_name)) {
          loaded_data[[file_name]] <- tryCatch(read.csv(con, row.names = 1), error = function(e) NULL)
        } else if (grepl("\\.tsv|\\.txt", file_name)) {
          loaded_data[[file_name]] <- tryCatch(read.delim(con, row.names = 1), error = function(e) NULL)
        } else {
          cat(sprintf("Unsupported .gz file format for %s. Skipped.\n", file_name))
        }
        close(con)

      } else if (grepl("\\.csv$", file_name)) {
        loaded_data[[file_name]] <- tryCatch(read.csv(file_path, row.names = 1), error = function(e) NULL)

      } else if (grepl("\\.tsv$|\\.txt$", file_name)) {
        loaded_data[[file_name]] <- tryCatch(read.delim(file_path, row.names = 1), error = function(e) NULL)

      } else {
        cat(sprintf("Unsupported format for %s. Skipping.\n", file_name))
      }
    }
  }

  if (read_to_memory) {
    varname <- paste0(gse_id, "_data")
    assign(varname, loaded_data, envir = .GlobalEnv)
    cat(sprintf("\nData loaded into .GlobalEnv as variable: %s\n", varname))
  } else {
    cat("\nDownload complete. Files saved only.\n")
  }

  invisible(NULL)
}
