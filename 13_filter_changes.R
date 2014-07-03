source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

library(rgdal)
library(stringr)
library(lubridate)
library(tools)

redo_classify <- FALSE
redo_extract <- FALSE
overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

sitecodes <- 'PSH'

image_basedir <- file.path(lcluc_folder, 'LCLUC_Classifications')

chgtraj_files <- dir(image_basedir,
                     pattern=paste0('^', sitecode, '_[0-9]{4}_[0-9]{4}_chgdetect_chgtraj.tif$'),
                     full.names=TRUE)
chgtraj_files <- dir(image_basedir,
                     pattern=paste0('^', sitecode, '_[0-9]{4}_[0-9]{4}_chgdetect_chgtraj.tif$'),
                     full.names=TRUE)

chgtraj_file <- chgtraj_files[1]
traj_freq <- foreach (chgtraj_file=iter(chgtraj_files),
                     .packages=c('teamlucc')) %do% {
    traj_freqs <- data.frame(freq(raster(chgtraj_file))
    chg_freqs <- chg_traj_out$lut
    chg_freqs$freq <- traj_freqs$count[match(chg_freqs$Code, traj_freqs$value)]
    chg_freqs <- chg_freqs[order(chg_freqs$t0_name),]
    freqs_filename <- file.path(image_basedir, paste(out_basename, 'chgtraj_freqs.csv', sep='_'))
    write.csv(chg_freqs, file=freqs_filename, row.names=FALSE)
}

}
