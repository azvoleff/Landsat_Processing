source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

library(rgdal)
library(stringr)
library(lubridate)
library(tools)

redo_chg_detection <- TRUE
overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

image_basedir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')
for (sitecode in sitecodes) {
    message(paste0('Performing change detection for ', sitecode, '...'))
    classes_files <- dir(image_basedir,
                         pattern=paste0('_predclasses.tif$'),
                         full.names=TRUE)
    probs_files <- dir(image_basedir,
                       pattern=paste0('_predprobs.tif$'),
                       full.names=TRUE)
    key_files <- dir(image_basedir,
                     pattern='_classeskey.csv$', full.names=TRUE)

    if (length(classes_files) == 0) {
        next
    }
    if (length(classes_files) == 1) {
        message(paste("cannot process", sitecode, "- only 1 classes file found"))
        next
    }

    classes_files_years <- gsub('_', '', str_extract(classes_files, 
                                                        '_[0-9]{4}_'))
    probs_files_years <- gsub('_', '', str_extract(probs_files, 
                                                      '_[0-9]{4}_'))
    key_files_years <- gsub('_', '', str_extract(key_files, '_[0-9]{4}_'))
    stopifnot(classes_files_years == probs_files_years)
    stopifnot(classes_files_years == key_files_years)
    file_years <- year(as.Date(key_files_years, '%Y'))

    # Check that all the key files are identical:
    code_keys <- lapply(key_files, read.csv)
    stopifnot(all(unlist(lapply(code_keys, identical, code_keys[[1]]))))
    class_key <- read.csv(key_files[1])
    classnames <- class_key$class

    # Setup pairs of image files
    classes_file_1s <- classes_files[1:(length(classes_files) - 1)]
    probs_file_1s <- probs_files[1:(length(probs_files) - 1)]
    probs_file_2s <- probs_files[2:length(probs_files)]
    year_1s <- file_years[1:(length(file_years) - 1)]
    year_2s <- file_years[2:length(file_years)]

    # ### TEMP ###
    # t1_classes <- stack(classes_file_1s)
    # t1_probs <- stack(probs_file_1s)
    # t2_probs <- stack(probs_file_2s)
    # year_1 <- year_1s[1]
    # year_2 <- year_2s[1]
    # output_path <- image_basedir
    # ### TEMP ###
    
    # Run change detection on each pair
    ret <- foreach (classes_file_1=iter(classes_file_1s), 
             probs_file_1=iter(probs_file_1s),
             probs_file_2=iter(probs_file_2s), 
             year_1=iter(year_1s), year_2=iter(year_2s), 
             .packages=c('teamlucc')) %do% {
        out_basename <- paste0(sitecode, '_', year_1, '-', year_2, 
                               '_chgdetect')

        output_files <- dir(image_basedir,
                            pattern=paste0('_chgtraj.tif$'),
                            full.names=TRUE)
        if (length(output_files) >= 1 & !redo_chg_detection) {
            next
        }

        t1_classes <- stack(classes_file_1)
        t1_probs <- stack(probs_file_1)
        t2_probs <- stack(probs_file_2)

        chg_dir_filename <- file.path(image_basedir,
                                      paste(out_basename, 'chgdir.tif', 
                                            sep='_'))
        chg_dir_image <- chg_dir(t1_probs, t2_probs, filename=chg_dir_filename, 
                                 overwrite=overwrite)

        chg_mag_filename <- file.path(image_basedir,
                                      paste(out_basename, 'chgmag.tif', 
                                            sep='_'))
        chg_mag_image <- chg_mag(t1_probs, t2_probs, filename=chg_mag_filename, 
                                 overwrite=overwrite)

        chg_threshold <- threshold(chg_mag_image, by=.025)

        chg_traj_filename <- file.path(image_basedir,
                                       paste(out_basename, 'chgtraj.tif', sep='_'))
        chg_traj_out <- chg_traj(t1_classes, chg_mag_image, chg_dir_image, 
                                 chg_threshold=chg_threshold,
                                 filename=chg_traj_filename,
                                 classnames=classnames,
                                 overwrite=overwrite)
    }
}
