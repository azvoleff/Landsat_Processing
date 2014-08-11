source('0_settings.R')

library(foreach)
library(itertools)
library(doParallel)

cl <- makeCluster(n_cpus)
registerDoParallel(cl)

library(rgdal)
library(raster)
library(ggplot2)
library(plyr)
library(grid)
library(stringr)
library(tools)

img_width <- 7
img_height <- 7
img_dpi <- 300

overwrite <- TRUE

zoi_folder <- file.path(prefix, 'TEAM', 'ZOIs')

image_basedir <- file.path(prefix, 'Landsat', 'Composites', 'Change_Detection')
out_dir <- image_basedir

stopifnot(file_test('-d', out_dir))

class_names_pretty <- c('Urban/built',
                        'Agriculture',
                        'Plantation forest',
                        'Natural forest',
                        'Other vegetation',
                        'Bare',
                        'Water',
                        'Unknown')
class_names_R <- c('Urban.built',
                   'Agriculture',
                   'Plantation.forest',
                   'Natural.forest',
                   'Other.vegetation',
                   'Bare',
                   'Water',
                   'Unknown')
class_names_abbrev <- c('Urban',
                        'Ag',
                        'PlanFor',
                        'NatFor',
                        'OthVeg',
                        'Bare',
                        'Water',
                        'Unk')

class_colors <- c('#CC0000',
                  '#F3F781',
                  '#3366FF',
                  '#088A08',
                  '#82FA58',
                  '#DBA901',
                  '#58D3F7',
                  '#A4A4A4')

chgtraj_files <- dir(image_basedir,
                     pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj.tif$'),
                     full.names=TRUE)
chgtraj_lut_files <- dir(image_basedir,
                         pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj_lut.csv$'),
                         full.names=TRUE)
stopifnot(length(chgtraj_files) == length(chgtraj_lut_files))

chgtraj_file <- chgtraj_files[1]
chgtraj_lut_file <- chgtraj_lut_files[1]

#' Function to plot classified image for a given year
#'
#' @export
#' @import rgdal
#' @import ggplot2
#' @importFrom plyr join
#' @importFrom grid unit
#' @importFrom sp spTransform CRS proj4string
#' @param x a forest change raster layer (a single layer of the layer 
#' stack output by \code{\link{annual_stack}}
#' @param aoi one or more AOI polygons as a \code{SpatialPolygonsDataFrame} 
#' object.  If there is a 'label' field  in the dataframe, it will be used to 
#' label the polygons in the plots. If the AOI is not in WGS 1984 (EPSG:4326), 
#' it will be reprojected to WGS84.
#' @param classes a \code{data.frame} with "code", "label", and (optionally) 
#' "color" columns. The "code" column indicates the numeric code in the image 
#' referring to a particular category or cover type, "label" indicates the 
#' label to use for each code, and "color" indicates the color to use on the 
#' image. "color" must be specified in a format recognized by \code{ggplot2}.
#' @param title_string the plot title
#' @param size_scale a number used to scale the size of the plot text
#' @param maxpixels the maximum number of pixels from x to use in plotting
plot_trajs <- function(x, aoi, classes, title_string='', size_scale=1, 
                         maxpixels=1e6, legend_title="Cover") {
    aoi_tr <- spTransform(aoi, CRS(proj4string(x)))
    aoi_tr$ID <- row.names(aoi_tr)
    if (!('label' %in% names(aoi_tr))) {
        aoi_tr$label <- paste('AOI', seq(1:nrow(aoi_tr@data)))
    }
    aoi_tr@data$id <- rownames(aoi_tr@data)
    aoi_points <- fortify(aoi_tr, region="id")
    aoi_df <- join(aoi_points, aoi_tr@data, by="id")

    if (ncell(x) > maxpixels) {
        x <- sampleRegular(x, maxpixels, asRaster=TRUE, useGDAL=TRUE)
    }
    dat <- as.data.frame(x, xy=TRUE)
    names(dat)[3] <- 'value'
    dat$value <- ordered(dat$value, levels=classes$code)

    if (!"color" %in% names(classes)) {
        gg_color_hue <- function(n) {
            # From: http://bit.ly/1yXZOlv
            hues = seq(15, 375, length=n+1)
            hcl(h=hues, l=65, c=100)[1:n]
        }
        classes$color <- gg_color_hue(nrow(classes))
    }

    long=lat=value=label=ID=NULL # For R CMD CHECK
    ggplot(dat) +
        geom_raster(aes(x, y, fill=value)) + coord_fixed() + 
        scale_fill_manual(legend_title, values=classes$color, breaks=classes$code,
                          labels=classes$label, drop=FALSE) +
        theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background=element_blank(), panel.border=element_blank(),
              panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              plot.background=element_blank(), axis.ticks=element_blank(),
              plot.margin=unit(c(.1, .1, .1, .1), 'cm')) +
        geom_path(aes(long, lat, group=ID), color="black", data=aoi_df, size=.7, alpha=.7) +
        ggtitle(title_string)
}

chgmag_file <- gsub("chgtraj", "chgmag", chgtraj_file)
chgmag_image <- raster(chgmag_file)
chgdir_file <- gsub("chgtraj", "chgdir", chgtraj_file)
chgdir_image <- raster(chgdir_file)
classes_1_probs_file <- "H:/Data/Landsat/Composites/Predictions/BIF_mosaic_1990_predictors_predprobs.tif"
classes_1_probs_image <- stack(classes_1_probs_file)
classes_2_probs_file <- "H:/Data/Landsat/Composites/Predictions/BIF_mosaic_1995_predictors_predprobs.tif"
classes_2_probs_image <- stack(classes_2_probs_file)

img1 <- as.matrix(c(.2, .4, .5), ncol=1, byrow=TRUE)
img2 <- as.matrix(c(.3, .1, .6), nrow=1, byrow=TRUE)
calc_chg_dir(as.matrix(c(.2, .4, .5), nrow=1, byrow=TRUE),
             as.matrix(c(.3, .1, .6), nrow=1, byrow=TRUE))

retvals <- foreach (chgtraj_file=iter(chgtraj_files), 
                    chgtraj_lut_file=iter(chgtraj_lut_files),
                    .packages=c('stringr', 'tools', 'raster', 'plyr', 'grid', 
                                'rgdal', 'ggplot2')) %dopar% {
    sitecode <- str_extract(basename(chgtraj_file), '^[a-zA-Z]*')
    year <- str_extract(chgtraj_file, '[0-9]{4}')
    chgtraj_rast <- raster(chgtraj_file)
    chgtraj_rast <- sampleRegular(chgtraj_rast, 1e6, asRaster=TRUE, 
                                  useGDAL=TRUE)

    # Mask out area outside of ZOI
    zoi_file <- dir(zoi_folder, pattern=paste0('^ZOI_', sitecode, '_[0-9]{4}.RData'), 
                    full.names=TRUE)
    stopifnot(length(zoi_file) == 1)
    load(zoi_file)
    zoi <- spTransform(zoi, CRS(proj4string(chgtraj_rast)))
    chgtraj_rast <- crop(chgtraj_rast, zoi)
    zoi_rast <- rasterize(zoi, chgtraj_rast, 1, silent=TRUE)
    chgtraj_rast[is.na(zoi_rast)] <- NA

    traj_codes <- read.csv(chgtraj_lut_file)
    traj_codes$t0_name_abbrev <- class_names_abbrev[match(traj_codes$t0_name, class_names_R)]
    traj_codes$t1_name_abbrev <- class_names_abbrev[match(traj_codes$t1_name, class_names_R)]
    traj_codes$trans_name <- with(traj_codes, paste(t0_name_abbrev, t1_name_abbrev, sep=' -> '))

    # Make forest/non-forest transition image
    traj_codes$ForestChange <- NA
    traj_codes$ForestChange[with(traj_codes, (t0_name != "Natural.forest") & (t1_name == "Natural.forest"))] <- 1
    traj_codes$ForestChange[with(traj_codes, (t0_name == "Natural.forest") & (t1_name != "Natural.forest"))] <- 2
    forest_subs <- data.frame(from=traj_codes$Code[!is.na(traj_codes$ForestChange)],
                              to=traj_codes$ForestChange[!is.na(traj_codes$ForestChange)])
    forest_change <- chgtraj_rast
    forest_change <- subs(chgtraj_rast, forest_subs)

    plot_trajs(forest_change, zoi,
               data.frame(code=c(1, 2), label=c("To forest", "From forest")),
               title="Forest transition")
    
    # Make loss/gain images
    trajs <- data.frame(code=NA,
                        label=class_names_abbrev,
                        color=class_colors,
                        stringsAsFactors=FALSE)

    plot_classes(chgtraj_rast_chgonly, zoi,
                 classes=data.frame(code=traj_codes$NewCode,
                                    label=traj_codes$trans_name),
                 legend_title="Transition")

    # Assign codes to trajs that are in this image. Assign codes to the 
    # others just so the plot code doesn't choke, and so these trajs remain 
    # in the legend, even if they do not appear in the image.
    #
    # First fill in proper codes for trajs that do appear in image
    trajs_rows_match <- match(traj_freqs$class, class_names_R)
    trajs[trajs_rows_match, ]$code <- traj_freqs$code
    # Now fill in random codes for trajs that don't appear in the image
    trajs_rows_nas <- is.na(trajs$code)
    trajs[trajs_rows_nas, ]$code <- seq(100, length.out=sum(trajs_rows_nas))

    title_string <- paste(sitecode, year, sep=' - ')
    plot_trajs(chgtraj_rast, zoi, trajs, title_string)
    ggsave(paste0(file_path_sans_ext(chgtraj_file), '_browse.png'), 
           width=img_width, height=img_height, dpi=img_dpi)
}

stopCluster(cl)
