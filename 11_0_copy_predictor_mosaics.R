source('0_settings.R')

library(stringr)

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

sitecodes <- c('YAS')

output_dir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')

stopifnot(file_test('-d', output_dir))

for (sitecode in sitecodes) {
    message(paste0('Copying files for ', sitecode, '...'))
    image_basedir <- file.path(prefix, 'Landsat', sitecode)
    orig_files <- dir(image_basedir,
                      pattern='^[a-zA-Z]*_mosaic_[0-9]{4}(_masks)?.tif',
                      full.names=TRUE)
    if (length(orig_files) == 0) {
        next
    }

    file_copies <- file.path(output_dir, basename(orig_files))

    file.copy(orig_files, file_copies, overwrite=TRUE)
}

for (sitecode in sitecodes) {
    message(paste0('Copying DEMs for ', sitecode, '...'))
    image_basedir <- file.path(prefix, 'Landsat', sitecode)
    orig_files <- dir(image_basedir,
                      pattern='^[a-zA-Z]*_mosaic_((dem)|(slopeaspect)).tif',
                      full.names=TRUE)
    if (length(orig_files) == 0) {
        next
    }

    file_copies <- file.path(output_dir, basename(orig_files))

    file.copy(orig_files, file_copies, overwrite=TRUE)
}
