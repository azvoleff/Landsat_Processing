source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

library(rgdal)
library(stringr)
library(lubridate)
library(tools)

redo_classify <- TRUE
overwrite <- TRUE

predictor_names <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7', 'msavi', 
                     'msavi_glcm_mean', 'msavi_glcm_variance', 
                     'msavi_glcm_dissimilarity', 'elev', 'slope', 'aspect')

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

image_basedir <- file.path(lcluc_folder, 'LCLUC_Classifications')
for (sitecode in sitecodes) {
    message(paste0('Processing ', sitecode, '...'))

    image_files <- dir(image_basedir,
                       pattern=paste0('^', sitecode, '_mosaic_[0-9]{4}_predictors.tif$'))

    output_files <- paste0(file_path_sans_ext(image_files), '_classified.tif')

    if (length(image_files) >= 1 & !overwrite) {
        image_files <- image_files[!file_test('-f', output_files)]
    }

    if (length(image_files) == 0) {
        next
    }

    model_file <- file.path(image_basedir, paste0(sitecode, '_rfmodel.RData'))
    if (!file_test('-f', model_file)) {
        next
    }
    load(model_file)

    ##########################################################################
    # Classify images
    for (image_file in image_files) {
        message(paste('Classifying', image_file))
        image_stack <- stack(file.path(image_basedir, image_file))
        fmask <- raster(file.path(image_basedir, 
                                 paste0(file_path_sans_ext(image_file), 
                                        '_masks', extension(image_file))), 
                        layer=2)
        # Assign standardized layer names to input image so that different 
        # images can be used with the same model
        names(image_stack) <- predictor_names

        results <- classify(image_stack, model)

        # Make a mask of 1s and NAs, with clear areas marked with 1s
        image_mask <- overlay(image_stack[[1]], fmask, fun=function(img, msk) {
            ret <- !is.na(img)
            ret[!ret] <- NA
            ret[(msk == 2) | (msk == 4) | (is.na(msk)) | (msk == 255)] <- NA
            return(ret)
        })
        
        classes <- results$classes * image_mask
        out_base <- file_path_sans_ext(file.path(image_basedir, image_file))
        classes_file <- paste0(out_base, '_predclasses', extension(image_file))
        classes <- writeRaster(classes, filename=classes_file, 
                               datatype='INT2S', overwrite=overwrite)

        probs <- round(results$probs * 100)
        probs <- probs * image_mask
        probs_file <- paste0(out_base, '_predprobs', extension(image_file))
        probs <- writeRaster(probs, filename=probs_file, datatype='INT2S', 
                             overwrite=overwrite)

        key_file <- paste0(out_base, '_classeskey.csv')
        write.csv(results$codes, file=key_file, row.names=FALSE)
    }
}

notify('Classification finished.')
