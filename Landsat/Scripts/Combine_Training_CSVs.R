combine_train_csvs <- function(
  input_dir = file.path("Input", "HighRes_Cover"),
  output_file = file.path("Output", "ShrubCover_R", "HighRes_Cover_Train_Combined.csv")
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }

  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }

  train_files <- sort(list.files(
    input_dir,
    pattern = "_Train\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  ))

  if (length(train_files) == 0L) {
    stop("No *_Train.csv files found under: ", input_dir)
  }

  column_names <- lapply(train_files, function(path) {
    names(data.table::fread(path, nrows = 0L, check.names = FALSE))
  })
  expected_names <- column_names[[1L]]
  schema_matches <- vapply(column_names, identical, logical(1), expected_names)

  if (!all(schema_matches)) {
    mismatched_files <- train_files[!schema_matches]
    stop(
      "CSV columns do not match the first file. Mismatched files: ",
      paste(mismatched_files, collapse = ", ")
    )
  }

  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  training_tables <- lapply(train_files, function(path) {
    data.table::fread(
      path,
      check.names = FALSE,
      showProgress = interactive()
    )
  })
  combined_data <- data.table::rbindlist(
    training_tables,
    use.names = TRUE,
    fill = FALSE
  )
  rm(training_tables)

  data.table::fwrite(combined_data, file = output_file)

  message("Combined ", length(train_files), " files and ", nrow(combined_data),
          " rows into: ", output_file)
  combined_data
}

if (sys.nframe() == 0L) {
  combined_train_data <- combine_train_csvs()
}
