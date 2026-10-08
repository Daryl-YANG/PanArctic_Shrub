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

terrain_folder <- 'Input/TerrainFeatures'
ccdc_folder <- "Input/CCDC_timeSeries"
save_folder <- "Input/ARD"
# create output directory if not exist
if (! file.exists(save_folder)) dir.create(save_folder,recursive=TRUE)
#*****************************************************************************************#


#*************************************** load data ***************************************#
#* load training database
topo_files <- list.files(terrain_folder, 
                         pattern = 'Stack_v2.tif$',
                         full.names = TRUE, 
                         recursive = TRUE)
#*****************************************************************************************#

# set up parallel backend (workers process CCDC files within each tile)
n_cores <- 4
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

for (f in topo_files)
{
  topo_rast <- terra::rast(f)
  
  filename <- basename(f)
  pattern <- gsub('_ArcticDEM_Terrain_Stack_v2.tif', '', filename)
  
  ccdc_loc <- paste0(ccdc_folder, '/', pattern)
  ccdc_files <- list.files(ccdc_loc,
                           pattern = '\\.tif$',
                           full.names = TRUE,
                           recursive = TRUE)

  if (length(ccdc_files) == 0) {
    warning('No CCDC rasters found for ', pattern, '; skipping.')
    next
  }
  
  ccdc_template <- terra::rast(ccdc_files[1])
  # Resample terrain onto the CCDC grid (CCDC CRS), then trim to its valid footprint
  topo_projected <- terra::project(topo_rast, ccdc_template, threads = TRUE)
  topo_projected <- terra::trim(topo_projected)

  # Spatial overlap of the two rasters, both now in the CCDC CRS
  overlap_ext <- terra::intersect(terra::ext(topo_projected), terra::ext(ccdc_template))
  if (is.null(overlap_ext)) {
    warning('No spatial overlap between terrain and CCDC for ', pattern, '; skipping.')
    next
  }
  topo_projected <- terra::crop(topo_projected, overlap_ext)
  topo_valid <- terra::app(topo_projected, function(values) {
    as.integer(all(!is.na(values)))
  })
  
  
  foldername <- paste0(save_folder, '/', pattern)
  if (! file.exists(foldername)) dir.create(foldername,recursive=TRUE)

  # SpatRasters cannot be sent to workers directly, so wrap them and unwrap per worker
  topo_wrapped <- terra::wrap(topo_projected)
  valid_wrapped <- terra::wrap(topo_valid)
  ext_vec <- as.vector(overlap_ext)

  status <- foreach(ccdc_file = ccdc_files, .packages = c("terra"),
                    .export = c("topo_wrapped", "valid_wrapped", "ext_vec", "foldername"),
                    .errorhandling = "pass") %dopar%
  {
    terra::terraOptions(threads = 1)
    topo_projected <- terra::unwrap(topo_wrapped)
    topo_valid <- terra::unwrap(valid_wrapped)
    overlap_ext <- terra::ext(ext_vec[["xmin"]], ext_vec[["xmax"]],
                              ext_vec[["ymin"]], ext_vec[["ymax"]])

    ccdc_rast <- terra::rast(ccdc_file)
    ccdc_rast <- try(terra::crop(ccdc_rast, overlap_ext), silent = TRUE)
    
    # Check if the result is an error; if so, skip this item
    if (inherits(ccdc_rast, "try-error")) {
      return(paste("Error encountered with item:", ccdc_file, "- Skipping!"))
    }

    comb_rast <- try(terra::mask(c(ccdc_rast/10000, topo_projected),
                 topo_valid,
                 maskvalues = 0), silent = TRUE)
    
    if (inherits(comb_rast, "try-error")) {
      return(paste("Error encountered with item:", ccdc_file, "- Skipping!"))
    }
    
    outname <- gsub('.tif', '_ARD.tif', basename(ccdc_file))
    outname <- paste0(foldername, '/', outname)
    terra::writeRaster(comb_rast, outname, overwrite = TRUE,
          gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9",
               "TILED=YES", "BIGTIFF=IF_SAFER"))
    # free worker memory before the next file
    rm(topo_projected, topo_valid, ccdc_rast, comb_rast)
    gc(verbose = FALSE)
    NULL
  }

  # report skipped or failed files
  for (i in seq_along(status)) {
    if (inherits(status[[i]], "error")) {
      cat("Error encountered with item:", ccdc_files[i], "-", conditionMessage(status[[i]]), "\n")
    } else if (is.character(status[[i]])) {
      cat(status[[i]], "\n")
    }
  }

  # free tile-level objects and temp files before the next tile
  rm(topo_rast, ccdc_template, topo_projected, topo_valid,
     topo_wrapped, valid_wrapped, ext_vec, status, ccdc_files)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE)
}

parallel::stopCluster(cl)
foreach::registerDoSEQ()
