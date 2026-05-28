

###########################################################################################
#
#        this script is for creating canopy height model from sfm or lidar point clouds
#
#    --- Last updated:  2023.02.15 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#


#********************************** support functions ************************************#
# function for extracting chm from point clouds 
chmExtr <- function (pcPATH, outDIR, reducer, classify)
{
  # create a temporary folder under the point could data path for parallel temporary outputs
  temp.dir <- paste0(pcPATH, '/', 'temp')
  if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)
  
  # load in point clouds (las or laz) file as catalog for parallel processing
  ctg <- readLAScatalog(data.dir, recursive = TRUE)
  
  ### set up parallel processing, please change these as necessary
  future::plan(multisession, workers = 6) 
  set_lidr_threads(6)
  opt_chunk_size(ctg) <- 100
  opt_chunk_buffer(ctg) <- 10
  
  ### split laz files into new tiles. this is not necessary if the original laz files are
  ### already tiles
  opt_output_files(ctg) <- paste0(temp.dir, "/retile_{XLEFT}_{YBOTTOM}")
  newctg = catalog_retile(ctg)
  
  ### reload the new tiles and check
  ctg <- readLAScatalog(temp.dir)
  las_check(ctg)
  
  
  ### reduce point cloud density if needed
  if (reducer == 'YES')
  {
    opt_output_files(ctg) <- paste0(temp.dir, "/{*}_thinned")
    thinned_ctg <- decimate_points(ctg, homogenize(dens, 1))
  } else 
  {
    thinned_ctg <- ctg
  }
  
  ### classify point cloud if needed
  if (classify == 'YES')
  {
    ### user can tweak the parameters as needed
    ws <- seq(3, 30, 6)
    th <- seq(0.1, 1.5, length.out = length(ws))
    opt_output_files(thinned_ctg) <- paste0(temp.dir, "/{*}_classified")
    classified_ctg <- classify_ground(thinned_ctg, algorithm = pmf(ws = ws, th = th))
  } else
  {
    # reassign classified ctg to ctg
    classified_ctg <- thinned_ctg
  }
  
  ### create digital terrain model raster ile
  opt_output_files(classified_ctg) <-  paste0(temp.dir, "/{*}_dtm")
  opt_stop_early(classified_ctg) <- FALSE
  dtm <- rasterize_terrain(classified_ctg, reso, tin())
  
  ### create canopy height laz file
  opt_output_files(classified_ctg) <-  paste0(temp.dir, "/{*}_norm")
  opt_stop_early(classified_ctg) <- FALSE
  ctg_norm <- normalize_height(classified_ctg, dtm)
  
  ### create canopy height raster file
  opt_output_files(ctg_norm) <- paste0(temp.dir, "/chm_{*}")
  opt_stop_early(ctg_norm) <- FALSE
  opt_merge(ctg_norm) <- TRUE
  chm_rst <- rasterize_canopy(ctg_norm, reso, p2r(0.15), overwrite=TRUE)
  
  ### remove NA values
  fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
  w <- matrix(1, 5, 5)
  chm_rst <- terra::focal(chm_rst, w, fun = fill.na)

  # calculate shrub cover
  opt_output_files(ctg_norm) <- paste0(temp.dir, "/shrubcover_{*}")
  #opt_output_files(ctg_norm) <- FALSE
  opt_stop_early(ctg_norm) <- FALSE
  opt_merge(ctg_norm) <- TRUE
  shrub_cover <- pixel_metrics(ctg_norm,~metrics(Z), res = 1, overwrite=TRUE)
  shrub_cover <- terra::rast(shrub_cover)
  
  # save intermediate files
  dtm_filename <- 'DTM.tif'
  dtm_filename <- paste0(outDIR, '/', dtm_filename)
  writeRaster(dtm, dtm_filename, overwrite = TRUE)
  
  chm_filename <- 'CHM.tif'
  chm_filename <- paste0(outDIR, '/', chm_filename)
  writeRaster(chm_rst, chm_filename, overwrite = TRUE)
  
  ### remove the temp folder
  unlink(temp.dir, recursive = T, force = T)
  # return result as a list
  return(list(chmPS = ctg_norm, chmRST = chm_rst, shrubCover = shrub_cover))
}

# function for calculating gcc from RGB imagery 
gccExtr <- function(rgbRST)
{
  # only extract the rgb bands
  rgb.rst <- rgbRST[[1:3]]
  names(rgb.rst) <- c('blue', 'green', 'red')
  # calculate gcc imagery
  gcc <- (rgb.rst$green)/app(rgb.rst[[1:3]], sum, cores = 6, na.rm = TRUE)
  gcc[rgb.rst$blue == 255 & rgb.rst$green == 255 & rgb.rst$red == 255] <- NA
  # return gcc
  return(gcc)
}

# define canopy cover function
metrics <- function(z)
{
  canopy=sum(z>0.4 & z<6) #anything over 0.4 meter
  perc=canopy/length(z)
  return(perc)
}
#*****************************************************************************************#


#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("lidR", 'future', 'terra')  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#


#************************************ user parameters ************************************#
# define input data directory
data.dir <- '/Users/dyd/Library/CloudStorage/Dropbox/Projects/NGEEArctic/Phase4/Codes/Test_Data'
### Create output folders
out.dir <- paste0(data.dir, '/', 'ShrubCover')
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)

### gcc related parameters
gccUSE <- "FALSE" #TRUE or FALSE
if(gccUSE == 'TRUE')
{
  gcc_thres = 0.4
}

### point cloud related parameters
# necessary when super dense cloud is not necessary, but not affect the result
reducer = 'YES'
if (reducer == 'YES')
{
  dens = 100 # point/m2
}

### define if you need to classify the point clouds. YES for classification, NO for don't
### classify (which means the point cloud is already classified)
classify = 'YES'

### define output resolution of terrain and canopy height raster
reso = 1
#*****************************************************************************************#

#************************************** main function ************************************#

### process point clouds data and return 
# chmPS: canopy height point cloud
# chmRST: canopy height rasterized file
# shrubCover: shrub cover derived directly from chmPS
ps_result <- chmExtr(data.dir, out.dir, reducer, classify)

# generate multi-scale shrub cover depending on the use of GCC or not
if(gccUSE == 'TRUE')
{
  ### calculate gcc from RBG imagery
  # load in RGB file
  
  rgb.dir <- list.files(data.dir, pattern = 'RGB.tif$',
                        recursive = TRUE, full.names = TRUE)
  rgb.rst <- terra::rast(rgb.dir)
  
  #extract gcc 
  gcc.rst <- gccExtr(rgb.rst)
  
  # save gcc file
  gcc.filename <- 'GCC.tif'
  gcc.filename <- paste0(out.dir, '/', gcc.filename)
  writeRaster(gcc.rst, gcc.filename, overwrite = TRUE)
  
  # remove rgb raster file from code to make space for point cloud processing
  rm(rgb.rst)
  
  ### convert gcc index into initial shrub cover for further processing
  gcc.rst[gcc.rst >= gcc_thres] <- 1
  gcc.rst[gcc.rst < gcc_thres] <- 0
  
  # calculate shrub canopy cover
  shrub.pix <- gcc.rst
  window.size <- round(reso/xres(shrub.pix))
  shrub.fcover <- aggregate(shrub.pix, 
                            fact = window.size, 
                            fun = mean, 
                            cores = 6,
                            na.rm = TRUE)
  
  rm(gcc.rst)
  # save chm files
  shrubcover.filename <- 'ShrubCover_Unfilterred.tif'
  shrubcover.filename <- paste0(out.dir, '/', shrubcover.filename)
  writeRaster(shrub.fcover, shrubcover.filename, overwrite = TRUE)
  
  # get chm raster from the point cloud processing result
  chm.rst <- ps_result$chmRST
  ### remove values that are considered as tree or any value below 0
  chm.rst[chm.rst < 0] <- 0
  chm.rst[chm.rst > 6] <- 0
  
  # clean GCC derived shrub cover based on canopy height
  chm.rst <- terra::resample(chm.rst, shrub.fcover)
  shrub.fcover[chm.rst < 0.4] <- 0
}
if(gccUSE == 'FALSE')
{  ### extract shrub cover derived from point could
  shrub.fcover <- ps_result$shrubCover
}

# resample shrub cover to different resolutions
shrubcover.filename <- 'ShrubCover_Clean_Orig.tif'
shrubcover.filename <- paste0(out.dir, '/', shrubcover.filename)
writeRaster(shrub.fcover, shrubcover.filename, overwrite = TRUE)

# resample shrub cover to 5 m resolution
scale_factor <- 5/xres(shrub.fcover)
shrub.fcover.5m <- aggregate(shrub.fcover, 
                             fact = scale_factor, 
                             fun = mean,
                             na.rm = TRUE)
shrubcover.filename <- 'ShrubCover_Clean_5m.tif'
shrubcover.filename <- paste0(out.dir, '/', shrubcover.filename)
writeRaster(shrub.fcover.5m, shrubcover.filename, overwrite = TRUE)

# resample shrub cover to 5 m resolution
scale_factor <- 30/xres(shrub.fcover)
shrub.fcover.30m <- aggregate(shrub.fcover, 
                              fact = scale_factor, 
                              fun = mean, 
                              na.rm = TRUE)
shrubcover.filename <- 'ShrubCover_Clean_30m.tif'
shrubcover.filename <- paste0(out.dir, '/', shrubcover.filename)
writeRaster(shrub.fcover.30m, shrubcover.filename, overwrite = TRUE)
#*****************************************************************************************#


#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#



