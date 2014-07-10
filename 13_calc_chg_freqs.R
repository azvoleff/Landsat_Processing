source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(4)

library(stringr)
library(tools)

redo_classify <- FALSE
redo_extract <- FALSE
overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

image_basedir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')

chgtraj_files <- dir(image_basedir,
                     pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj.tif$'),
                     full.names=TRUE)
chgtraj_lut_files <- dir(image_basedir,
                         pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj_lut.csv$'),
                         full.names=TRUE)
stopifnot(length(chgtraj_files) == length(chgtraj_lut_files))

traj_freq <- foreach (chgtraj_file=iter(chgtraj_files),
                      chgtraj_lut_file=iter(chgtraj_lut_files),
                      .packages=c('teamlucc', 'tools')) %do% {
    traj_freqs <- data.frame(freq(raster(chgtraj_file)))
    chg_freqs <- read.csv(chgtraj_lut_file)
    chg_freqs$freq <- traj_freqs$count[match(chg_freqs$Code, traj_freqs$value)]
    chg_freqs <- chg_freqs[order(chg_freqs$t0_name, chg_freqs$t1_name),]
    freqs_filename <- paste0(file_path_sans_ext(chgtraj_file), '_freqs.csv')
    write.csv(chg_freqs, file=freqs_filename, row.names=FALSE)
}
