source('0_settings.R')

library(stringr)
library(tools)
library(dplyr)
library(reshape2)

library(ggplot2)

library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(n_cpus)

overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

zoi_folder <- file.path(prefix, 'TEAM', 'ZOIs')
image_basedir <- file.path(prefix, 'Landsat', 'LCLUC_Classifications')

reprocess <- TRUE
imgtype <- 'raw'

###############################################################################
# Summarize mosaics
mosaic_stats_RData_file <- 'mosaic_pixel_stats.RData'
if (reprocess) {
    stopifnot(imgtype %in% c('normalized', 'raw'))
    if (imgtype == 'normalized') {
        pattern <- '^[a-zA-Z]*_mosaic_normalized_[0-9]{4}.tif$'
    } else {
        pattern <- '^[a-zA-Z]*_mosaic_[0-9]{4}.tif$'
    }
    mosaic_files <- dir(image_basedir, pattern=pattern, full.names=TRUE, recursive=TRUE)

    sitecodes <- str_extract(basename(mosaic_files), '^[a-zA-Z]*')
    mosaic_files <- mosaic_files[sitecodes != 'BBS']

    mosaic_stats <- foreach(mosaic_file=iter(mosaic_files),
                            .packages=c('raster', 'tools', 'stringr'),
                            .combine=rbind) %dopar% {
        fmask <- raster(paste0(file_path_sans_ext(mosaic_file), '_masks', 
                                   extension(mosaic_file)),
                        band=2)
        sitecode <- str_extract(basename(mosaic_file), '^[a-zA-Z]*')
        year <- as.numeric(str_extract(basename(mosaic_file), '[0-9]{4}'))

        # Mask out area outside of ZOI
        zoi_file <- dir(zoi_folder, pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                        full.names=TRUE)
        stopifnot(length(zoi_file) == 1)
        load(zoi_file)
        zoi <- spTransform(zoi, CRS(proj4string(fmask)))

        # Set the update value to 99 so that NA values inside the mask can be 
        # properly picked up, and NAs outside the mask can be ignored (since 
        # they are coded with the ignored 99 code.
        fmask <- mask(fmask, zoi, value=99)

        bs <- blockSize(fmask)
        num_fill <- 0
        num_cloud <- 0
        num_clear <- 0
        num_water <- 0
        num_snow <- 0
        num_missing <- 0
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
        num_pixels <- num_fill + num_cloud + num_clear + num_water + num_snow + num_missing
        return(data.frame(site=sitecode, date=year, num_fill=num_fill, 
                          num_missing=num_missing, num_cloud=num_cloud, 
                          num_clear=num_clear, num_water=num_water, 
                          num_snow=num_snow, total_pixels=num_pixels))
    }

    # Add rows with num_fill equal to the number of pixels in the whole image for 
    # sites without any data for a particular year.  Do this with a merge.
    miss_data <- data.frame(site=rep(unique(mosaic_stats$site), each=5))
    miss_data$date <- rep(unique(mosaic_stats$date), length.out=nrow(miss_data))
    miss_data$num_fill <- mosaic_stats$total_pixels[match(miss_data$site, 
                                                          mosaic_stats$site)]
    miss_data$total_pixels <- miss_data$num_fill
    miss_data <- miss_data[!(paste(miss_data$site, miss_data$date) %in% 
                             paste(mosaic_stats$site, mosaic_stats$date)), ]

    stopifnot(sum(is.na(mosaic_stats)) == 0)
    mosaic_stats <- merge(mosaic_stats, miss_data, all=TRUE)
    mosaic_stats[is.na(mosaic_stats)] <- 0
    save(mosaic_stats, file=mosaic_stats_RData_file)
}
load(mosaic_stats_RData_file)

missing_summary <- summarize(group_by(mosaic_stats, site, date),
                             pct_cloud=(num_cloud / (num_clear + num_water + num_cloud + num_snow))*100,
                             pct_missing=(num_fill / (num_fill + num_clear + num_water + num_cloud + num_snow)*100))
#NaN results from mosaicks with ALL data missing
missing_summary$pct_cloud[is.nan(missing_summary$pct_cloud)] <- 0
filter(missing_summary, pct_cloud > 1)
filter(missing_summary, pct_cloud > 5)
filter(missing_summary, pct_missing > 5)

plot_codes <- data.frame(site=unique(missing_summary$site))
plot_codes$shape <- factor(rep(1:4, length.out=nrow(plot_codes)))
plot_codes$colour <- factor(rep(1:4, each=4, length.out=nrow(plot_codes)))
plot_codes$linetype <- factor(rep(1:4, each=4, length.out=nrow(plot_codes)))
mosaic_stats <- merge(mosaic_stats, plot_codes)

missing_summary <- merge(missing_summary, plot_codes)

ggplot(missing_summary) +
    geom_line(aes(date, pct_cloud, colour=colour, linetype=linetype, group=site)) +
    geom_point(aes(date, pct_cloud, colour=colour, shape=shape, group=site)) +
    scale_colour_manual("Site", labels=plot_codes$site, breaks=plot_codes$colour, values=mosaic_stats$colour) +
    scale_shape_manual("Site", labels=plot_codes$site, breaks=plot_codes$shape, values=mosaic_stats$shape) +
    scale_linetype_manual("Site", labels=plot_codes$site, breaks=plot_codes$linetype, values=mosaic_stats$linetype)

ggplot(missing_summary) +
    geom_line(aes(date, pct_missing, colour=colour, linetype=linetype, group=site)) +
    geom_point(aes(date, pct_missing, colour=colour, shape=shape, group=site)) +
    scale_colour_manual("Site", labels=plot_codes$site, breaks=plot_codes$colour, values=mosaic_stats$colour) +
    scale_shape_manual("Site", labels=plot_codes$site, breaks=plot_codes$shape, values=mosaic_stats$shape) +
    scale_linetype_manual("Site", labels=plot_codes$site, breaks=plot_codes$linetype, values=mosaic_stats$linetype)



cloud_wide_table <- dcast(missing_summary, site ~ date)
cloud_wide_table[2:ncol(cloud_wide_table)] <- round(cloud_wide_table[2:ncol(cloud_wide_table)], 2)
write.csv(cloud_wide_table, file='mosaic_pixel_missing_summary.csv', row.names=FALSE)

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
        these_files <- dir(image_basedir, pattern=pattern)
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
