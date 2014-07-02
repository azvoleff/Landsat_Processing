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

for (sitecode in sitecodes) {
    message(paste0('Filtering classification for ', sitecode, '...'))

    image_basedir <- file.path(prefix, 'Landsat', sitecode)
    image_files <- dir(image_basedir,
                       pattern='^[a-zA-Z]*_mosaic_[0-9]{4}_predictors.tif$')



    raster

    # Run change detection on each pair
    ret <- foreach (chgtraj_file=iter(chgtraj_files),
                    .packages=c('teamlucc')) %do% {
        traj_freqs <- data.frame(freq(chg_traj_out$traj))
        chg_freqs <- chg_traj_out$lut
        chg_freqs$freq <- traj_freqs$count[match(chg_freqs$Code, traj_freqs$value)]
        chg_freqs <- chg_freqs[order(chg_freqs$t0_name),]
        freqs_filename <- file.path(image_basedir, paste(out_basename, 'chgtraj_freqs.csv', sep='_'))
        write.csv(chg_freqs, file=freqs_filename, row.names=FALSE)
    }

}
