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
overwrite <- TRUE

predictor_names <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b7', 'msavi', 
                     'msavi_glcm_mean', 'msavi_glcm_variance', 
                     'msavi_glcm_dissimilarity', 'elev', 'slope', 'aspect')

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

image_files <- c()
model_files <- c()
image_basedir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')
for (sitecode in sitecodes) {
    these_image_files <- dir(image_basedir,
                       pattern=paste0('^', sitecode, '_mosaic_[0-9]{4}_predictors.tif$'))

    output_files <- paste0(file_path_sans_ext(these_image_files),
                           '_predclasses', extension(image_file))

    if (length(these_image_files) >= 1 & !redo_classify) {
        these_image_files <- these_image_files[!file_test('-f', output_files)]
    }

    if (length(these_image_files) == 0) {
        next
    }

    this_model_file <- file.path(image_basedir, paste0(sitecode, '_rfmodel.RData'))
    if (!file_test('-f', this_model_file)) {
        next
    }

    image_files <- c(image_files, these_image_files)
    model_files <- c(model_files, rep(this_model_file, length(these_image_files)))
}

stopifnot(length(image_files) == length(model_files))

notify(paste0('Starting classification. ', length(image_files), ' images to process.'))
num_res <- foreach (image_file=iter(image_files),
                    model_file=iter(model_files),
                    .packages=c('teamlucc', 'tools', 'stringr', 'notifyR'),
                    .inorder=FALSE) %dopar% {
    raster_tmpdir <- paste0(tempdir(), '_raster_',
                            paste(sample(c(letters, 0:9), 15), collapse=''))
    dir.create(raster_tmpdir)
    rasterOptions(tmpdir=raster_tmpdir)
    
    sitecode <- str_extract(basename(image_file), '^[a-zA-Z]')

    load(model_file)

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

    removeTmpFiles(h=0)
    unlink(raster_tmpdir)

    notify(paste0('Finished classifying ', image_file))
    return(1)
}

if (length(num_res) == 0) num_res <- 0
notify(paste0('Classification finished. Classified ', sum(num_res), ' images.'))
