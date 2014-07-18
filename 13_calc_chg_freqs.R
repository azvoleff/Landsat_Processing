source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

library(stringr)
library(tools)

redo_classify <- FALSE
redo_extract <- FALSE
overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

zoi_folder <- file.path(prefix, 'TEAM', 'ZOIs')
image_basedir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')

chgtraj_files <- dir(image_basedir,
                     pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj.tif$'),
                     full.names=TRUE)
chgtraj_lut_files <- dir(image_basedir,
                         pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj_lut.csv$'),
                         full.names=TRUE)
stopifnot(length(chgtraj_files) == length(chgtraj_lut_files))

foreach (chgtraj_file=iter(chgtraj_files),
         chgtraj_lut_file=iter(chgtraj_lut_files),
         .packages=c('tools', 'stringr', 'raster', 'rgdal')) %dopar% {
    sitecode <- str_extract(basename(chgtraj_file), '^[a-zA-Z]*')

    chgtraj_rast <- raster(chgtraj_file)

    # Mask out area outside of ZOI
    zoi_file <- dir(zoi_folder, pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                    full.names=TRUE)
    stopifnot(length(zoi_file) == 1)
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(proj4string(chgtraj_rast)))

    # Set masked areas to 255 (code for fill)
    chgtraj_rast <- mask(chgtraj_rast, zoi, update=255)

    traj_freqs <- data.frame(freq(chgtraj_rast))
    chg_freqs <- read.csv(chgtraj_lut_file)
    chg_freqs$freq <- traj_freqs$count[match(chg_freqs$Code, traj_freqs$value)]
    chg_freqs <- chg_freqs[order(chg_freqs$t0_name, chg_freqs$t1_name),]
    freqs_filename <- paste0(file_path_sans_ext(chgtraj_file), '_freqs.csv')
    write.csv(chg_freqs, file=freqs_filename, row.names=FALSE)
}

# Output simple class frequency tables for each image
predclasses_files <- dir(image_basedir,
                         pattern=paste0('^[a-zA-Z]{2,3}_mosaic_[0-9]{4}_predictors_predclasses.tif$'),
                         full.names=TRUE)
classeskey_files <- dir(image_basedir,
                        pattern=paste0('^[a-zA-Z]{2,3}_mosaic_[0-9]{4}_predictors_classeskey.csv$'),
                        full.names=TRUE)
stopifnot(length(classeskey_files) == length(predclasses_files))
foreach (predclasses_file=iter(predclasses_files),
         classeskey_file=iter(classeskey_files),
         .packages=c('stringr', 'tools', 'raster')) %dopar% {
    sitecode <- str_extract(basename(predclasses_file), '^[a-zA-Z]*')
    year <- str_extract(predclasses_file, '[0-9]{4}')
    classes_rast <- stack(predclasses_file)

    # Mask out area outside of ZOI
    zoi_file <- dir(zoi_folder, pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                    full.names=TRUE)
    stopifnot(length(zoi_file) == 1)
    load(zoi_file)

    # Set masked areas to 255 (code for fill)
    classes_rast <- mask(classes_rast, zoi, update=255)

    classeskey <- read.csv(classeskey_file)
    class_freqs <- data.frame(freq(classes_rast))
    names(class_freqs) <- c('code', 'freq')
    class_freqs <- cbind(sitecode=sitecode, year=year,
                         name=classeskey$class[match(class_freqs$code, classeskey$code)],
                         class_freqs)
    write.csv(class_freqs, file=paste0(file_path_sans_ext(predclasses_file), 
                                       '_classfreqs.csv'))
}
