source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

cl <- makeCluster(n_cpus)
registerDoParallel(cl)

library(rgdal)
library(stringr)
library(lubridate)
library(tools)

redo_chg_detection <- TRUE
overwrite <- TRUE

min_chg_threshold <- 20

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code
sitecodes <- c("BIF", "CAX", "COU", "CSN",
               "MAS", "PSH", "RNF", "VB",
               "YAN", "YAS", "BCI", "BBS",
               "UDZ", "NAK")

zoi_folder <- file.path(prefix, 'TEAM', 'ZOIs')
predictions_dir <- file.path(prefix, 'Landsat', 'Composites', 'Predictions')
chgdetect_dir <- file.path(prefix, 'Landsat', 'Composites', 'Change_Detection')
out_dir <- file.path(prefix, 'Landsat', 'Composites', 'Change_Detection')

chgmag_files <- c()
zoi_files <- c()
for (sitecode in sitecodes) {
    these_chgmag_files <- dir(chgdetect_dir,
                               pattern=paste0('^', sitecode, '_[0-9]{4}-[0-9]{4}_chgdetect_chgmag.tif$'),
                               full.names=TRUE)

    if (length(these_chgmag_files) == 0) {
        next
    }
    this_zoi_file <- dir(zoi_folder,
                         pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                         full.names=TRUE)
    stopifnot(length(this_zoi_file) == 1)
    
    chgmag_files <- c(chgmag_files, these_chgmag_files)
    
    zoi_files <- c(zoi_files, rep(this_zoi_file, length(these_chgmag_files)))
}
stopifnot(length(chgmag_files) == length(zoi_files))

# Run change detection on each pair
notify(paste0('Starting change detection. ', length(chgmag_files), ' image sets to process.'))
num_res <- foreach (chgmag_file=iter(chgmag_files), zoi_file=iter(zoi_files),
                    .packages=c('teamlucc', 'notifyR', 'stringr', 'rgdal'),
                    .combine=c, .inorder=FALSE) %dopar% {
    raster_tmpdir <- file.path(temp, paste0('raster_',
                               paste(sample(c(letters, 0:9), 15), collapse='')))
    dir.create(raster_tmpdir)
    rasterOptions(tmpdir=raster_tmpdir)

    sitecode <- str_extract(basename(chgmag_file), '^[a-zA-Z]*')

    year_1 <- as.numeric(gsub('[_-]','', str_extract(chgmag_file, '_[0-9]{4}-')))
    year_2 <- as.numeric(gsub('[_-]','', str_extract(chgmag_file, '-[0-9]{4}_')))
    stopifnot(year_1 < year_2)

    out_basename <- paste0(sitecode, '_', year_1, '-', year_2, '_chgdetect')

    output_files <- dir(out_dir, pattern=paste0(out_basename, '_chgtraj.tif'), 
                        full.names=TRUE)
    if (length(output_files) >= 1 & !redo_chg_detection) {
        return()
    }

    classes_1_filename <- dir(predictions_dir,
                              pattern=paste0(sitecode, '_mosaic_', year_1,
                                             '_predictors_predclasses.tif'),
                              full.names=TRUE)
    stopifnot(length(classes_1_filename) == 1)
    classes_1_image <- stack(classes_1_filename)

    chgmag_image <- raster(chgmag_file)

    chgdir_filename <- file.path(out_dir,
                                  paste(out_basename, 'chgdir.tif', 
                                        sep='_'))
    chgdir_image <- raster(chgdir_filename)

    # Load ZOI for use in masking images
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(proj4string(chgmag_image)))
    zoi_rast <- rasterize(zoi, chgmag_image, 1, silent=TRUE)


    chgmag_image_crop <- crop(chgmag_image, zoi)
    chgmag_image_crop <- mask(chgmag_image_crop, zoi)
    chg_threshold <- threshold(chgmag_image_crop)
    
    chg_threshold_filename <- file.path(out_dir,
                                        paste(out_basename, 'threshold.txt', 
                                              sep='_'))
    write.csv(data.frame(sitecode=sitecode, year_1=year_1, year_2=year_2, 
                         threshold=chg_threshold), file=chg_threshold_filename, 
              row.names=FALSE)

    if (chg_threshold < min_chg_threshold) chg_threshold <- min_chg_threshold
    
    key_file_1 <- gsub('predclasses.tif', 'classeskey.csv', classes_1_filename)
    class_key <- read.csv(key_file_1)
    classnames <- class_key$class

    chg_traj_out <- chg_traj(chgmag_image, chgdir_image, 
                             chg_threshold=chg_threshold,
                             classnames=classnames)

    chg_traj_filename <- file.path(out_dir, paste(out_basename, 'chgtraj.tif', 
                                                  sep='_'))
    # Recheck - but below masking line should no longer be needed
    #chg_traj_image <- chg_traj_out$traj * image_mask
    chg_traj_image <- writeRaster(chg_traj_out, 
                                  filename=chg_traj_filename, 
                                  overwrite=overwrite, datatype='INT2S')

    chg_traj_lut_filename <- file.path(out_dir,
                                       paste(out_basename, 'chgtraj_lut.csv', 
                                             sep='_'))
    lut <- traj_lut(class_key$code, class_key$class)
    write.csv(lut, file=chg_traj_lut_filename, row.names=FALSE)

    # Set masked areas to 99 so they can be differentiated. Don't use mask as 
    # it has a bug where it doesn't set NA areas in the image to the 
    # updatevalue
    chg_traj_image_masked <- chg_traj_image
    chg_traj_image_masked[is.na(zoi_rast)] <- 99

    traj_freqs <- data.frame(freq(chg_traj_image_masked))
    chg_freqs <- lut
    chg_freqs$freq <- traj_freqs$count[match(chg_freqs$Code, traj_freqs$value)]
    chg_freqs <- chg_freqs[order(chg_freqs$t0_name, chg_freqs$t1_name),]
    freqs_filename <- file.path(out_dir, paste(out_basename, 
                                               'chgtraj_freqs.csv', sep='_'))
    write.csv(chg_freqs, file=freqs_filename, row.names=FALSE)

    removeTmpFiles(h=0)
    unlink(raster_tmpdir)

    return(1)
}

if (length(num_res) == 0) num_res <- 0
notify(paste0('Finished change detection. Processed ', sum(num_res), ' images.'))

stopCluster(cl)
