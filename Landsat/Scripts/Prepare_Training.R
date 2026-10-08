###########################################################################################
#
#  this script extracts spectral reflectance (and quality control layer if preferred) from
#  original modis brdf-corrected reflectance hdf files downloaded from DACC, convert hdf 
#  to tif and merge spectral bands. the output is a single raster file for each hdf file
#  that contains all desired spectral bands and quality control layers
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
list.of.packages <- c("ggplot2", 'terra', 'stringr', 'sf')  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
#* define work dir
setwd("/Volumes/Utility_SSD/NGEEArctic/panArctic_Shrub")

shrubcover_folder <- 'Input/HighRes_Cover'
ccdc_folder <- 'Input/CCDC_timeSeries'
topo_folder <- 'Input/TerrainFeatures'

output_folder <- 'Output/ShrubCover_R'
#* create output directory if not exist
if (! file.exists(output_folder)) dir.create(output_folder,recursive=TRUE)
#*****************************************************************************************#

#*********************************** load files ******************************************#
#* search all refereence shrub cover files from shrubcover_folder
shrubcover_files <- list.files(shrubcover_folder,
                               pattern = '*ShrubCover_30m.tif$',
                               recursive = TRUE,
                               full.names = TRUE)

#* print out a preview of the files
str(shrubcover_files)
#* summarize how many files were found
print(paste0('no. of files found: ', length(shrubcover_files)))


ccdc_files <- list.files(ccdc_folder,
                         pattern = '*.tif$',
                         recursive = TRUE,
                         full.names = TRUE)
#* print out a preview of the files
str(ccdc_files)
#* summarize how many files were found
print(paste0('no. of files found: ', length(ccdc_files)))


topo_files <- list.files(topo_folder,
                         pattern = '*.tif$',
                         recursive = TRUE,
                         full.names = TRUE)
#* print out a preview of the files
str(topo_files)
#* summarize how many files were found
print(paste0('no. of files found: ', length(topo_files)))

# Cache raster footprints so overlap checks do not reopen every candidate raster.
cache_footprints <- function(files) {
  setNames(lapply(files, function(f) {
    sf::st_as_sfc(sf::st_bbox(terra::rast(f)))
  }), files)
}

ccdc_footprints <- cache_footprints(ccdc_files)
topo_footprints <- cache_footprints(topo_files)

find_overlapping_files <- function(files, footprints, shrubcover_box) {
  files[vapply(files, function(f) {
    candidate_box <- footprints[[f]]
    shrubcover_box_in_crs <- sf::st_transform(
      shrubcover_box,
      sf::st_crs(candidate_box)
    )
    any(sf::st_intersects(shrubcover_box_in_crs, candidate_box,
                          sparse = FALSE))
  }, logical(1))]
}
#*****************************************************************************************#
# 1. Expand memory allowance if working with large files
terraOptions(memfrac = 0.8, memmax = 80) # e.g., allow up to 12 GB
# 3. Increase GDAL cache size (e.g., 64GB)
gdalCache(size = 65536)

#*********************************** main function ***************************************#
### go through all reference shrub cover files and find overlapping CCDC reference image 
#* and topo image, then extract predictor variables
for (filename in shrubcover_files)
{
  print(filename)
  #* identify the year of reference data collection
  yr_coll <- str_extract(filename, "\\d{4}")
  
  #* load in reference shrub cover raster layer
  shrubcover_rast <- terra::rast(filename)
  names(shrubcover_rast) <- 'shrubcover'
  shrubcover_box <- sf::st_as_sfc(sf::st_bbox(shrubcover_rast))
  
  ###* filter ccdc files by collection year before checking spatial overlap
  str_search <- paste0('v20240207_', yr_coll)
  year_ccdc_files <- ccdc_files[grepl(str_search, ccdc_files, fixed = TRUE)]
  tp_overlapping_ccdc_files <- find_overlapping_files(
    year_ccdc_files, ccdc_footprints, shrubcover_box
  )

  ###* search for topo files that are spatially overlap with reference
  sp_overlapping_topo_files <- find_overlapping_files(
    topo_files, topo_footprints, shrubcover_box
  )
  
  
  if (length(sp_overlapping_topo_files) >= 1)
  {
    raster_list <- lapply(sp_overlapping_topo_files, rast)
    raster_collection <- sprc(raster_list)
    topo_rast <- mosaic(raster_collection, fun = "mean")
  } else (next)
  
  if (length(tp_overlapping_ccdc_files) >= 1)
  {
    raster_list <- lapply(tp_overlapping_ccdc_files, rast)
    raster_collection <- sprc(raster_list)
    ccdc_rast <- mosaic(raster_collection, fun = "mean")
  } else (next)
  
  
  # Reproject or match CRS if necessary
  if (crs(ccdc_rast) != crs(shrubcover_rast)) {
    shrubcover_rast <- project(shrubcover_rast, crs(ccdc_rast),
                               threads = TRUE, use_gdal = TRUE)
  }
  #* clip ccdc files
  ccdc_clip <- terra::crop(ccdc_rast, shrubcover_rast)
  
  
  #topo_rast <- terra::rast(sp_overlapping_topo_files)
  # Reproject or match CRS if necessary
  if (crs(ccdc_rast) != crs(topo_rast)) {
    topo_rast <- project(topo_rast, crs(ccdc_rast),
                         threads = TRUE, use_gdal = TRUE)
  }
  #* clip topo files
  topo_clip <- terra::crop(topo_rast, shrubcover_rast)
  
  
  #* make sure all raster files on the same resolution and extent
  shrubcover_rast <- resample(shrubcover_rast, ccdc_clip)
  topo_clip <- resample(topo_clip, ccdc_clip)
  
  merged_rast <- c(shrubcover_rast, ccdc_clip/10000, topo_clip)
  
  df <- as.data.frame(merged_rast, na.rm = TRUE)
  
  file_basename <- basename(filename)
  df <- data.frame('site' = rep(file_basename, nrow(df)),
                     df)

  outname <- gsub('ShrubCover_30m.tif', 'ShrubCover_30m_Train.csv', filename)
  write.csv(df, outname,row.names = FALSE)
  
  rm(ccdc_rast, topo_rast, topo_rast)
  rm(df)
}
































