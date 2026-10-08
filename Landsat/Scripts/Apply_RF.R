###########################################################################################
#

#
#    --- Last updated:  2021.06.22 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#


#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", 
                      "ranger", 
                      "ggpmisc", 
                      "terra", 
                      "foreach", 
                      "doParallel")
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"]) ]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
# clean any lingering parallel backend from previous work in the same R session
if (exists(".foreachGlobals", mode = "list")) rm(".foreachGlobals", envir = .GlobalEnv)
if (exists("cl", mode = "list")) { try(stopCluster(cl), silent = TRUE); rm("cl", envir = .GlobalEnv) }
foreach::registerDoSEQ()
#*****************************************************************************************#



#************************************ user parameters ************************************#
#* define work dir
setwd("/Volumes/NGEE/NGEEArctic/panArctic_Shrub")

ard_folder <- 'Input/ARD'
models_folder <- "Output/ShrubCover_Map/models"
save_folder <- "Output/ShrubCover_Map/Maps"
# create output directory if not exist
if (! file.exists(save_folder)) dir.create(save_folder,recursive=TRUE)
#*****************************************************************************************#


#*************************************** load data ***************************************#
#* load training database
rf_models <- list.files(models_folder, 
                        pattern = '^fcover_model_.*shrubcover\\.rds$',
                        full.names = TRUE, 
                        recursive = TRUE)

if (length(rf_models) == 0L) {
  stop('No saved shrub-cover models found in ', models_folder, call. = FALSE)
}

ard_files <- list.files(ard_folder, 
                        pattern = 'ARD.tif$',
                        full.names = TRUE, 
                        recursive = TRUE)
# ignore macOS resource-fork files (._*) on external volumes
ard_files <- ard_files[!startsWith(basename(ard_files), '._')]

if (length(ard_files) == 0L) {
  stop('No ARD rasters found in ', ard_folder, call. = FALSE)
}
#*****************************************************************************************#
#************************************* apply model ***************************************#
# Cells per processing block; bounds memory use in each worker
chunk_cells <- 4e6

# One persistent worker per model. Each worker loads its model once (dropping training-only
# components) and keeps it for the whole run; ranger scales poorly across threads, so
# single-threaded workers running in parallel are much faster than one multi-threaded predict.
cluster <- parallel::makeCluster(length(rf_models))
invisible(parallel::clusterEvalQ(cluster, {
  library(terra)
  library(ranger)
  terra::terraOptions(threads = 1)
}))
worker_features <- parallel::clusterApply(cluster, rf_models, function(f) {
  obj <- readRDS(f)
  obj$model$inbag.counts <- NULL
  obj$model$predictions <- NULL
  assign('worker_obj', obj, envir = .GlobalEnv)
  obj$metadata$feature_names
})
feature_names <- worker_features[[1]]
stopifnot(all(vapply(worker_features, identical, logical(1), feature_names)))
# positions of the reflectance bands within the predictors
refl_idx <- grep('^spectral_(blue|green|red|NIR|SWIR1|SWIR2)_', feature_names)

# Runs on each worker: read one block of rows, keep only valid pixels (some reflectance band
# non-zero and no missing predictors), and predict them with the worker's model.
predict_block <- function(ard_file, feature_names, refl_idx, start, n) {
  r <- terra::rast(ard_file)[[feature_names]]
  terra::readStart(r)
  x <- terra::readValues(r, row = start, nrows = n, mat = TRUE)
  terra::readStop(r)
  valid <- rowSums(x[, refl_idx, drop = FALSE] != 0, na.rm = TRUE) > 0 &
    rowSums(is.na(x)) == 0L
  if (!any(valid)) {
    rm(r, x, valid)
    return(list(idx = integer(0), pred = numeric(0)))
  }
  obj <- get('worker_obj', envir = .GlobalEnv)
  idx <- (start - 1L) * terra::ncol(r) + which(valid)
  x <- x[valid, , drop = FALSE]
  pred <- predict(obj$model, data = x, num.threads = 1)$predictions
  # release this block's data in the worker; only the model stays resident
  rm(r, x, valid, obj)
  list(idx = idx, pred = pred)
}

# Free cached memory in every worker and in the main session so nothing carries over to the next tile.
clean_memory <- function(cluster) {
  invisible(parallel::clusterEvalQ(cluster, {
    gc(verbose = FALSE, full = TRUE)
    NULL
  }))
  invisible(gc(verbose = FALSE, full = TRUE))
}

# Ensemble mean and between-model SD for every pixel of one tile; invalid pixels stay NA.
predict_ensemble <- function(ard_file, cluster, feature_names, refl_idx, chunk_cells) {
  template <- terra::rast(ard_file)[[feature_names]]
  n_col <- terra::ncol(template)
  n_row <- terra::nrow(template)
  n_mod <- length(cluster)
  cover_mean <- rep(NA_real_, terra::ncell(template))
  cover_sd   <- rep(NA_real_, terra::ncell(template))
  rows_per_chunk <- max(1L, floor(chunk_cells / n_col))
  
  for (start in seq(1L, n_row, by = rows_per_chunk)) {
    n <- min(rows_per_chunk, n_row - start + 1L)
    res <- parallel::clusterCall(cluster, predict_block, ard_file, feature_names,
                                 refl_idx, start, n)
    idx <- res[[1]]$idx
    if (length(idx) == 0L) next
    pred_sum   <- Reduce(`+`, lapply(res, function(z) z$pred))
    pred_sumsq <- Reduce(`+`, lapply(res, function(z) z$pred^2))
    cover_mean[idx] <- pred_sum / n_mod
    cover_sd[idx] <- if (n_mod > 1L) {
      sqrt(pmax((pred_sumsq - pred_sum^2 / n_mod) / (n_mod - 1L), 0))
    } else NA_real_
  }
  list(mean = cover_mean, sd = cover_sd, template = template)
}

for (ard_file in ard_files)
{
  ard_id   <- tools::file_path_sans_ext(basename(ard_file))
  out_file <- file.path(save_folder, paste0(ard_id, '_shrubcover.tif'))
  # skip tiles that were already completed so interrupted runs can resume
  if (file.exists(out_file)) {
    message('Skipping (exists): ', ard_file)
    next
  }
  
  pred <- predict_ensemble(ard_file, cluster, feature_names, refl_idx, chunk_cells)
  
  cover_out <- terra::rast(pred$template, nlyrs = 2)
  terra::values(cover_out) <- cbind(pred$mean, pred$sd)
  names(cover_out) <- c('shrubcover_mean', 'shrubcover_sd')
  
  # write to a partial file first so an interrupted run never leaves a finished-looking tile
  partial_file <- file.path(save_folder, paste0(ard_id, '_shrubcover_partial.tif'))
  terra::writeRaster(cover_out, partial_file, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"))
  file.rename(partial_file, out_file)
  
  message('Finished ', ard_file)
  rm(pred, cover_out)
  clean_memory(cluster)
}
stopCluster(cluster)
