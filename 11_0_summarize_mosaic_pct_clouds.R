source('0_settings.R')

library(stringr)
library(tools)
library(dplyr)
library(reshape2)

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

reprocess <- FALSE

###############################################################################
# Summarize mosaics

mosaic_stats_RData_file <- 'mosaic_pixel_stats.RData'
if (reprocess) {
    imgtype <- 'raw'
    stopifnot(imgtype %in% c('normalized', 'raw'))
    base_dir <- file.path(prefix, 'Landsat')
    if (imgtype == 'normalized') {
        pattern <- '^[a-zA-Z]*_mosaic_normalized_[0-9]{4}.tif$'
    } else {
        pattern <- '^[a-zA-Z]*_mosaic_[0-9]{4}.tif$'
    }
    mosaic_files <- dir(base_dir, pattern=pattern, full.names=TRUE, recursive=TRUE)

    mosaic_stats <- foreach(mosaic_file=iter(mosaic_files),
                            .packages=c('raster', 'tools', 'stringr'),
                            .combine=rbind) %dopar% {
        fmask <- raster(paste0(file_path_sans_ext(mosaic_file), '_masks', 
                                   extension(mosaic_file)),
                        band=2)
        bs <- blockSize(fmask)
        num_missing <- 0
        num_fill <- 0
        num_cloud <- 0
        num_clear <- 0
        num_water <- 0
        num_snow <- 0
        for (block_num in 1:bs$n) {
            fmask_bl <- getValuesBlock(fmask, row=bs$row[block_num], 
                                       nrows=bs$nrows[block_num])
            num_missing <- num_missing + sum(is.na(fmask_bl))
            num_fill <- num_fill + sum(fmask_bl == 255, na.rm=TRUE)
            num_cloud <- num_cloud + sum(fmask_bl == 2, na.rm=TRUE) + sum(fmask_bl == 4, na.rm=TRUE)
            num_clear <- num_clear + sum(fmask_bl == 0, na.rm=TRUE)
            num_water <- num_water + sum(fmask_bl == 1, na.rm=TRUE)
            num_snow <- num_snow + sum(fmask_bl == 3, na.rm=TRUE)
        }
        num_pixels <- num_fill + num_missing + num_cloud + num_clear + num_water + num_snow
        sitecode <- str_extract(basename(mosaic_file), '^[a-zA-Z]*')
        year <- str_extract(basename(mosaic_file), '[0-9]{4}')
        return(data.frame(site=sitecode, date=year, num_fill=num_fill, 
                   num_missing=num_missing, num_cloud=num_cloud, 
                   num_clear=num_clear, num_water=num_water, num_snow=num_snow, 
                   total_pixels=num_pixels))
    }
    save(mosaic_stats, file=mosaic_stats_RData_file)
}

load(mosaic_stats_RData_file)

cloud_pcts <- summarize(group_by(mosaic_stats, site, date),
                        pct_cloud=(num_cloud / total_pixels)*100)
filter(cloud_pcts, pct_cloud > 1)

cloud_wide_table <- dcast(cloud_pcts, site ~ date)
cloud_wide_table[2:ncol(cloud_wide_table)] <- round(cloud_wide_table[2:ncol(cloud_wide_table)], 2)
write.csv(cloud_wide_table, file='mosaic_pixel_cloud_pcts.csv', row.names=FALSE)


###############################################################################
# Summarize DEM mosaics
#
# - Min elevation
# - Max elevation
# - Mean elevation
# - Mean slope
# - Max slope
get_dem_status <- function(pattern, name) {
    statuses <- foreach(sitecode=iter(sitecodes), .combine=rbind) %do% {
        base_dir <- file.path(prefix, 'Landsat', sitecode)
        these_files <- dir(base_dir, pattern=pattern)
        if (length(these_files) >= 1) {
            statuses <- data.frame(site=sitecode, this_status=TRUE)
        } else {
            statuses <- data.frame()
        }
        return(statuses)
    }
    names(statuses)[names(statuses) == 'this_status'] <- name
    return(statuses)
}

mosaic_dem_status <- get_dem_status('^[a-zA-Z]*_mosaic_dem.tif$', 'mosaic_dem')
mosaic_slopeaspect_status <- get_dem_status('^[a-zA-Z]*_mosaic_slopeaspect.tif$', 'mosaic_slpasp')
mosaic_dem_statuses <- merge(mosaic_dem_status, mosaic_slopeaspect_status, all=TRUE)
mosaic_dem_statuses 
