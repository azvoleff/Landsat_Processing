source('0_settings.R')

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

library(gdalUtils)
library(rgdal)
library(rgeos)
library(stringr)
library(tools)
library(lubridate)

reprocess <- FALSE
overwrite <- TRUE
builddem <- TRUE
#imgtype <- 'normalized'
imgtype <- 'raw'

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

sitecodes <- sitecodes[sitecodes != 'BBS']

stopifnot(imgtype %in% c('normalized', 'raw'))

for (sitecode in sitecodes) {
    message(paste0('Mosaicking images for ', sitecode, '...'))

    if (imgtype == 'normalized') {
        pattern='^[a-zA-Z]*_[0-9]{3}-[0-9]{3}_[0-9]{4}-[0-9]{3}_cf(_((normbase)|(normalized))).tif$'
    } else {
        pattern='^[a-zA-Z]*_[0-9]{3}-[0-9]{3}_[0-9]{4}-[0-9]{3}_cf.tif$'
    }
    base_dir <- file.path(prefix, 'Landsat', sitecode)
    image_files <- dir(base_dir, pattern=pattern, full.names=TRUE)
    image_stacks <- lapply(image_files, stack)

    mask_files <- paste0(file_path_sans_ext(image_files), '_masks.tif')
    if (any(!file_test('-f', mask_files))) {
        stop('could not locate mask files')
    }
    mask_stacks <- lapply(mask_files, stack)

    image_date_strings <- unique(str_extract(basename(image_files), 
                                             '_[0-9]{4}-[0-9]{3}_'))

    # Function to get maximal extent from a list of extents
    get_total_extent <- function(extents) {
        full_ext <- extents[[1]]
        extents <- extents[-1]
        while (length(extents) > 0) {
            full_ext <- merge(full_ext, extents[[1]])
            extents <- extents[-1]
        }
        return(full_ext)
    }
    # First calculate maximal extent of each stack of images, considering all 
    # the dates to ensure that the mosaic extents match for each epoch even if 
    # some epochs are missing path/rows
    mos_exts <- foreach(image_date_string=iter(image_date_strings),
                        .packages=c('raster', 'rgdal', 'lubridate', 'tools', 
                                    'foreach', 'iterators')) %dopar% {
        image_date_object <- as.Date(image_date_string, '_%Y-%j_')
        epoch_image_files <- image_files[grepl(image_date_string, image_files)]
        # Calculate full extent of mosaic
        extents <- lapply(epoch_image_files, function(x) {
            x <- stack(x)
            extent(x)
        })
        get_total_extent(extents)
    }
    mos_ext <- get_total_extent(mos_exts)

    # Now perform mosaicking
    mosaic_stacks <- foreach(image_date_string=iter(image_date_strings),
                             .packages=c('raster', 'rgdal', 'lubridate', 
                                         'tools', 'foreach', 'iterators',
                                         'gdalUtils'),
                             .combine=c) %dopar% {
        image_date_object <- as.Date(image_date_string, '_%Y-%j_')

        epoch_image_files <- image_files[grepl(image_date_string, image_files)]
        epoch_mask_files <- paste0(file_path_sans_ext(epoch_image_files), '_masks.tif')

        stopifnot(all(file_test('-f', epoch_image_files)))
        stopifnot(all(file_test('-f', epoch_mask_files)))

        epoch_year <- year(as.Date(image_date_string, '_%Y-%j_'))
        if (imgtype == 'normalized') {
            mosaic_out_file <- file.path(base_dir,
                                         paste0(sitecode, '_mosaic_normalized_', 
                                                year(image_date_object),  
                                            extension(image_files[1])))
        } else {
            mosaic_out_file <- file.path(base_dir,
                                         paste0(sitecode, '_mosaic_', 
                                                year(image_date_object),  
                                            extension(image_files[1])))
        }
        if (file_test('-f', mosaic_out_file) & !reprocess) {
            return()
        }

        mask_out_file <- paste0(file_path_sans_ext(mosaic_out_file), '_masks', 
                                extension(mosaic_out_file))

        epoch_images <- lapply(epoch_image_files, brick)
        epoch_masks <- lapply(epoch_mask_files, brick)

        masked_epoch_image_files <- foreach(epoch_image=iter(epoch_images), 
                                            epoch_mask=iter(epoch_masks), 
                                            .packages=c('raster', 'rgdal', 
                                                        'gdalUtils'),
                                            .combine=c) %do% {
            # Make sure clouds are masked out (NA) and make sure missing values 
            # and SLC-off is masked out.
            #
            # Remember layer 2 is fmask
            out_file <- extension(rasterTmpFile(), '.tif')
            epoch_image <- overlay(epoch_image, epoch_mask[[2]], 
                                   fun=function(img, msk) {
                img[(msk == 2) | (msk == 4) | (msk == 255)] <- NA 
                img[is.na(msk)] <- NA 
                return(img)
            }, datatype=dataType(epoch_image)[1], filename=out_file)
            return(out_file)
        }

        # Setup the images to have their origins at 0,0
        mosaic_te <- as.numeric(bbox(mos_ext))
        tr <- c(30, 30)
        # Setup xmin
        mosaic_te[1] <- round(mosaic_te[1] - mosaic_te[1] %% tr[1])
        # Setup ymin
        mosaic_te[2] <- round(mosaic_te[2] - mosaic_te[2] %% tr[2])
        # Setup xmax
        mosaic_te[3] <- round(mosaic_te[3] + tr[1] - mosaic_te[3] %% tr[1])
        # Setup ymax
        mosaic_te[4] <- round(mosaic_te[4] + tr[2] - mosaic_te[4] %% tr[2])
        stopifnot(all(round(mosaic_te / 30) == (mosaic_te / 30)))

        mask_stack <- gdalwarp(epoch_mask_files,
                               dstfile=mask_out_file,
                               r='near', output_Raster=TRUE, 
                               of='GTiff', dstnodata="None",
                               overwrite=overwrite, multi=TRUE, 
                               wo=paste0("NUM_THREADS=", n_cpus), 
                               te=mosaic_te, tr=c(30, 30),
                               ot='Byte')

        image_stack <- gdalwarp(masked_epoch_image_files,
                                dstfile=mosaic_out_file,
                                r='cubicspline', output_Raster=TRUE, 
                                of='GTiff',
                                overwrite=overwrite, multi=TRUE, 
                                wo=paste0("NUM_THREADS=", n_cpus), 
                                te=mosaic_te, tr=c(30, 30),
                                ot='Int16')
    }

    # Check extents of all mosaics are equal
    if (imgtype == 'normalized') {
        pattern <- '^[a-zA-Z]*_mosaic_normalized_[0-9]{4}.tif$'
    } else {
        pattern <- '^[a-zA-Z]*_mosaic_[0-9]{4}.tif$'
    }
    mosaic_files <- dir(base_dir, pattern=pattern, full.names=TRUE)
    mosaic_stacks <- lapply(mosaic_files, stack)
    mos_exts <- lapply(mosaic_stacks, extent)
    for (mos_ext in mos_exts) {
        stopifnot(mos_ext == mos_exts[[1]])
    }

    dem_mosaic_filename <- file.path(base_dir,
                                     paste0(sitecode, '_mosaic_dem.tif'))
    if (builddem & (!file_test('-f', dem_mosaic_filename) | reprocess)) {
        message(paste0('Mosaicking DEMs for ', sitecode, '...'))
        
        mos_ext <- as(mos_exts[[1]], 'SpatialPolygons')

        proj4string(mos_ext) <- proj4string(mosaic_stacks[[1]])

        mos_ext_dem_proj <- spTransform(mos_ext, CRS(proj4string(dem_extents)))

        intersecting <- as.logical(gIntersects(dem_extents, 
                                               gUnaryUnion(mos_ext_dem_proj), byid=TRUE))
        if (sum(intersecting) == 0) {
            stop('no intersecting dem extents found')
        }

        dem_list <- dem_extents[intersecting, ]$filename
        dem_rasts <- lapply(dem_list, raster)

        to_srs <- proj4string(mosaic_stacks[[1]])

        # Calculate minimum bounding box coordinates:
        dem_te <- as.numeric(bbox(mos_ext))
        to_res <- c(30, 30)
        dem_mosaic <- gdalwarp(dem_list, dstfile=dem_mosaic_filename,
                               te=dem_te, t_srs=to_srs, tr=to_res, 
                               r='cubicspline', output_Raster=TRUE, multi=TRUE, 
                               of='GTiff',
                               wo=paste0("NUM_THREADS=", n_cpus), 
                               overwrite=overwrite)

        # Note that the default output of 'terrain' is in radians
        slopeaspect <- terrain(dem_mosaic, opt=c('slope', 'aspect'))
        slopeaspect$aspect <- calc(slopeaspect$aspect, fun=function(vals) {
            vals[vals >= 2*pi] <- 0
            vals
            })
        # Note that slopeaspect is scaled - slope by 10000, and aspect by 1000 so 
        # that the layers can be saved as INT2S
        slopeaspect <- stack(round(raster(slopeaspect, layer=1) * 10000),
                             round(raster(slopeaspect, layer=2) * 1000))

        slopeaspect_mosaic_filename <- file.path(base_dir,
                                         paste0(sitecode, '_mosaic_slopeaspect.tif'))
        slopeaspect <- writeRaster(slopeaspect, filename=slopeaspect_mosaic_filename, 
                                 overwrite=overwrite, datatype='INT2S')
    }

}
