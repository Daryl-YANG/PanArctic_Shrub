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


#************************************ user functions *************************************#
# Sample an equal number of rows from each non-empty target bin.
equal_bin_sample <- function(df, target_col = 'yVar', max_n = 150000L, n_bins = 25L, seed = 42L) {
  if (nrow(df) <= max_n) {
    return(df)
  }
  
  if (max_n < 1L || n_bins < 1L) {
    stop("max_n and n_bins must both be positive.")
  }
  
  #set.seed(NULL)
  valid_rows <- is.finite(df[[target_col]])
  df <- df[valid_rows, , drop = FALSE]
  if (nrow(df) == 0) {
    return(df)
  }
  
  y <- df[[target_col]]
  y_min <- min(y)
  y_max <- max(y)
  if (y_max <= y_min) {
    return(df[sample.int(nrow(df), min(max_n, nrow(df))), , drop = FALSE])
  }
  
  breaks <- seq(y_min, y_max, length.out = n_bins + 1L)
  bin_id <- cut(y, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  bin_counts <- tabulate(bin_id, nbins = n_bins)
  populated_bins <- which(bin_counts > 0L)
  if (length(populated_bins) == 0L) {
    return(df[0, , drop = FALSE])
  }
  
  samples_per_bin <- min(
    floor(max_n / length(populated_bins)),
    min(bin_counts[populated_bins])
  )
  if (samples_per_bin < 1L) {
    stop("max_n is too small to draw at least one sample from each populated bin.")
  }
  
  sampled_idx <- unlist(lapply(populated_bins, function(bin) {
    bin_rows <- which(bin_id == bin)
    bin_rows[sample.int(length(bin_rows), samples_per_bin)]
  }), use.names = FALSE)
  df[sampled_idx, , drop = FALSE]
}

# JAB-style uncertainty from tree predictions: each tree is a bootstrap replicate,
# and the between-tree spread provides an uncertainty estimate for the prediction.
get_jackknife_uncertainty <- function(model, validation_data) {
  pred_all <- predict(model, data = validation_data, predict.all = TRUE)$predictions
  if (is.null(dim(pred_all))) {
    return(rep(0, nrow(validation_data)))
  }
  pred_mean <- rowMeans(pred_all, na.rm = TRUE)
  pred_sd <- apply(pred_all, 1, sd, na.rm = TRUE)
  list(mean = pred_mean, sd = pred_sd)
}
#*****************************************************************************************#


#************************************ user parameters ************************************#
#* define work dir
setwd("/Volumes/NGEE/NGEEArctic/panArctic_Shrub")

training_data_folder <- 'Input/Training_Data'
out_folder <- "Output/ShrubCover_Map"
# create output directory if not exist
if (! file.exists(out_folder)) dir.create(out_folder,recursive=TRUE)
# creat an temporary to store files temporarily generated during the course of processing
model_folder <- file.path(paste0(out_folder, "/", 'models'))
if (! file.exists(model_folder)) dir.create(model_folder,recursive=TRUE)

# define how many models you want to build?
nmodel <- 10

# large-data training settings: keep the distribution of the target while capping runtime
max_training_n <- 5000000L
train_fraction <- 0.8
hist_bins <- 25L

# pft names
yname <- 'shrubcover'
#*****************************************************************************************#


#*************************************** load data ***************************************#
#* load training database
train_database_dir <- list.files(training_data_folder, 
                                 pattern = 'HighRes_Cover_Train_Combined.csv$',
                                 full.names = TRUE, 
                                 recursive = TRUE)
data.org <- read.csv(train_database_dir)
str(data.org)
data.org <- na.omit(data.org)

#* split the database randomly into a calibration data and validation subsets
split_ind <- sample(seq_len(nrow(data.org)), 
                     size = nrow(data.org) * train_fraction)
cal_data <- data.org[split_ind, ]
val_data   <- data.org[-split_ind, ]
#*****************************************************************************************#


#************************************** train model **************************************#
### train random forest model to predict shrub cover from landsat reflectance and topo features
#* construct a number of models defined in nmodel to account for training data uncertainty
#* use ranger instead of RandomForest to increase efficiency
#* use doParallel cluster to perform paralell process

# prepare data for random forest 
xVAR <- cal_data[, 3:29, drop = FALSE]
rf.data <- data.frame(yVar = cal_data[[yname]], xVAR)
rf.data <- rf.data[complete.cases(rf.data), , drop = FALSE]

#Setup backend to use many processors
totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(10) #totalCores[1]-1
registerDoParallel(cluster)
result <- foreach(i=1:nmodel) %dopar% {
  data.train <- equal_bin_sample(rf.data, target_col = 'yVar', max_n = max_training_n,
                                           n_bins = hist_bins, seed = 42)
  fcover.model <- ranger::ranger(
    formula = yVar ~ ., 
    data = data.train,
    num.trees = 100,
    mtry = max(2L, floor(sqrt(ncol(data.train) - 1L))),
    min.node.size = 5,
    importance = "impurity",
    keep.inbag = TRUE,
    write.forest = TRUE,
    num.threads = 10
  )

  #* save the model
  model.out <- file.path(model_folder, paste0('fcover_model_',i, '_', yname, '.rds'))
  model.meta <- list(
    yname = yname,
    feature_names = names(xVAR),
    max_training_n = max_training_n,
    train_fraction = train_fraction,
    hist_bins = hist_bins,
    model_type = 'ranger_random_forest'
  )
  saveRDS(list(model = fcover.model, metadata = model.meta), model.out)
  
  pred_val <- predict(fcover.model, data = val_data)$predictions
}
stopCluster(cluster)

pred.df <- data.frame(result)

### calculate mean cover prediction
pred.mean <- as.data.frame(apply(pred.df, 1, FUN = mean, na.rm = TRUE))
names(pred.mean) <- 'pred_mean'
### calculate cover prediction uncertainty
pred.unc <- as.data.frame(apply(pred.df, 1, FUN = sd, na.rm = TRUE))
names(pred.unc) <- 'pred_unc'


plotdata <- data.frame('site' = val_data$site, 
                       'truth' = val_data$shrubcover,
                       'predicted' = pred.mean$pred_mean,
                       'unc' = pred.unc$pred_unc)

# calculate absolute prediction error
mae <- sum(abs(plotdata$predicted-plotdata$truth)/nrow(plotdata))
mae <- format(round(mae, 5), nsmall = 5)

formula <- y ~ x
ggplot(data = plotdata, aes(x = truth, y = predicted)) +
  geom_hex(bins = 100) +
  scale_fill_continuous(type = 'viridis', limits=c(1, 20), 
                        na.value = "yellow") +
  geom_smooth(method = 'lm', formula = formula, col = 'black') +
  geom_abline(intercept = 0, slope = 1, 
              linetype = 'dashed', 
              color = 'red',
              size = 1.5) +
  ylim(c(0, 1)) + xlim(c(0, 1)) +
  labs(x = 'Truth Cover', y = 'Predicted Cover') +
  theme(legend.position = 'none') +
  theme(axis.text = element_text(size=12), axis.title=element_text(size=13)) +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x = 0.65, label.y = 0.25,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x = 0.65, label.y = 0.2,
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  annotate('text', x= 0.66, y = 0.12, label = paste0('MAE = ', mae), size = 4, hjust = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


png.name <- paste0(out_folder, "/",'model_eva_hex.pdf')
ggsave(png.name, plot = last_plot(), width = 12, height = 12, units = 'cm')

















