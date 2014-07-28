source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(6)

library(rgdal)
library(stringr)
library(lubridate)
library(tools)

redo_chg_detection <- TRUE
overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code
sitecodes <- c('MAS')

zoi_folder <- file.path(prefix, 'TEAM', 'ZOIs')
image_basedir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')
classes_file_1s <- c()
classes_file_2s <- c()
zoi_files <- c()
for (sitecode in sitecodes) {
    these_classes_files <- dir(image_basedir,
                               pattern=paste0('^', sitecode, 
                                              '_mosaic_[0-9]{4}_predictors_predclasses.tif$'),
                               full.names=TRUE)

    if (length(these_classes_files) == 0) {
        next
    }
    if (length(these_classes_files) == 1) {
        message(paste("cannot process", sitecode, "- only 1 classes file found"))
        next
    }

    this_zoi_file <- dir(zoi_folder,
                         pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                         full.names=TRUE)
    stopifnot(length(this_zoi_file) == 1)
    
    # Setup pairs of image files
    these_classes_file_1s <- these_classes_files[1:(length(these_classes_files) - 1)]
    these_classes_file_2s <- these_classes_files[2:length(these_classes_files)]

    classes_file_1s <- c(classes_file_1s, these_classes_file_1s)
    classes_file_2s <- c(classes_file_2s, these_classes_file_2s)
    
    zoi_files <- c(zoi_files, rep(this_zoi_file, length(these_classes_files) - 1))
}
stopifnot(length(classes_file_1s) == length(classes_file_2s))
stopifnot(length(classes_file_1s) == length(zoi_files))

# Run change detection on each pair
notify(paste0('Starting change detection. ',
              length(classes_file_1s), ' image sets to process.'))
num_res <- foreach (classes_file_1=iter(classes_file_1s), 
                    classes_file_2=iter(classes_file_2s),
                    zoi_file=iter(zoi_files),
                    .packages=c('teamlucc', 'notifyR', 'stringr', 'rgdal'),
                    .combine=c, .inorder=FALSE) %dopar% {
    raster_tmpdir <- file.path(temp, paste0('raster_',
                               paste(sample(c(letters, 0:9), 15), collapse='')))
    dir.create(raster_tmpdir)
    rasterOptions(tmpdir=raster_tmpdir)

    sitecode <- str_extract(basename(classes_file_1), '^[a-zA-Z]*')

    year_1 <- as.numeric(gsub('_','', str_extract(classes_file_1,
                                                  '_[0-9]{4}_')))
    year_2 <- as.numeric(gsub('_','', str_extract(classes_file_2,
                                                  '_[0-9]{4}_')))
    stopifnot(year_1 < year_2)

    out_basename <- paste0(sitecode, '_', year_1, '-', year_2, 
                           '_chgdetect')

    output_files <- dir(image_basedir,
                        pattern=paste0('_chgtraj.tif$'),
                        full.names=TRUE)
    if (length(output_files) >= 1 & !redo_chg_detection) {
        return()
    }

    mask_file_1 <- gsub('predclasses', 'masks', classes_file_1)
    mask_file_2 <- gsub('predclasses', 'masks', classes_file_2)
    probs_file_1 <- gsub('predclasses', 'predprobs', classes_file_1)
    probs_file_2 <- gsub('predclasses', 'predprobs', classes_file_2)

    key_file_1 <- gsub('predclasses.tif', 'classeskey.csv', classes_file_1)
    key_file_2 <- gsub('predclasses.tif', 'classeskey.csv', classes_file_2)
    stopifnot(read.csv(key_file_1) == read.csv(key_file_2))
    class_key <- read.csv(key_file_1)
    classnames <- class_key$class

    # Make a mask of 1s and NAs, with clear areas marked with 1s. This needs to 
    # take into account BOTH images.
    mask_1 <- raster(mask_file_1, layer=2)
    mask_2 <- raster(mask_file_2, layer=2)
    image_mask <- overlay(mask_1, mask_2, fun=function(msk1, msk2) {
        ret <- !is.na(msk1)
        ret[(msk1 == 2) | (msk1 == 4) | (is.na(msk1)) | (msk1 == 255)] <- NA
        ret[(msk2 == 2) | (msk2 == 4) | (is.na(msk2)) | (msk2 == 255)] <- NA
        return(ret)
    })

    t1_classes <- stack(classes_file_1)
    t1_probs <- stack(probs_file_1)
    t2_probs <- stack(probs_file_2)

    chg_dir_filename <- file.path(image_basedir,
                                  paste(out_basename, 'chgdir.tif', 
                                        sep='_'))
    chg_dir_image <- chg_dir(t1_probs, t2_probs)
    chg_dir_image <- chg_dir_image * image_mask
    writeRaster(chg_dir_image, filename=chg_dir_filename, 
                overwrite=overwrite, datatype=dataType(chg_dir_image))

    chg_mag_filename <- file.path(image_basedir,
                                  paste(out_basename, 'chgmag.tif', 
                                        sep='_'))
    chg_mag_image <- chg_mag(t1_probs, t2_probs)
    chg_mag_image <- chg_mag_image * image_mask
    writeRaster(chg_mag_image, filename=chg_mag_filename, 
                overwrite=overwrite, datatype=dataType(chg_mag_image))
   
    chg_threshold <- threshold(chg_mag_image, by=.025)
    
    chg_traj_out <- chg_traj(t1_classes, chg_mag_image, chg_dir_image, 
                             chg_threshold=chg_threshold,
                             classnames=classnames)
    
    chg_traj_filename <- file.path(image_basedir,
                                   paste(out_basename, 'chgtraj.tif', sep='_'))
    chg_traj_image <- chg_traj_out$traj * image_mask
    chg_traj_image <- writeRaster(chg_traj_image, filename=chg_traj_filename, 
                                  overwrite=overwrite, datatype='INT2S')

    chg_traj_lut_filename <- file.path(image_basedir,
                                       paste(out_basename, 'chgtraj_lut.csv', 
                                             sep='_'))
    write.csv(chg_traj_out$lut, file=chg_traj_lut_filename, row.names=FALSE)

    # Calculate change frequencies, masking out area outside of ZOI
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(proj4string(chg_traj_image)))
    zoi <- rasterize(zoi, chg_traj_image, 1, silent=TRUE)

    # Set masked areas to 99 so they can be differentiated. Don't use mask as 
    # it has a bug where it doesn't set NA areas in the image to the 
    # updatevalue
    chg_traj_image_masked <- chg_traj_image
    chg_traj_image_masked[is.na(zoi)] <- 99

    traj_freqs <- data.frame(freq(chg_traj_image_masked))
    chg_freqs <- chg_traj_out$lut
    chg_freqs$freq <- traj_freqs$count[match(chg_freqs$Code, traj_freqs$value)]
    chg_freqs <- chg_freqs[order(chg_freqs$t0_name, chg_freqs$t1_name),]
    freqs_filename <- file.path(image_basedir,
                                paste(out_basename, 'freqs.csv', sep='_'))
    write.csv(chg_freqs, file=freqs_filename, row.names=FALSE)

    removeTmpFiles(h=0)
    unlink(raster_tmpdir)

    return(1)
}

if (length(num_res) == 0) num_res <- 0
notify(paste0('Finished change detection. Processed ', sum(num_res), ' images.'))
