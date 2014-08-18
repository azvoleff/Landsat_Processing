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
type <- "cairo"

overwrite <- TRUE

zoi_folder <- file.path(prefix, 'TEAM', 'ZOIs')

image_basedir <- file.path(prefix, 'Landsat', 'Composites', 'Change_Detection')
out_dir <- image_basedir

stopifnot(file_test('-d', out_dir))

chgtraj_files <- dir(image_basedir,
                     pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj.tif$'),
                     full.names=TRUE)
chgtraj_lut_files <- dir(image_basedir,
                         pattern=paste0('^[a-zA-Z]{2,3}_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj_lut.csv$'),
                         full.names=TRUE)
stopifnot(length(chgtraj_files) == length(chgtraj_lut_files))

# chgtraj_file <- chgtraj_files[1]
# chgtraj_lut_file <- chgtraj_lut_files[1]

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
                       maxpixels=2e6, legend_title="Cover") {
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
    ForestChange <- chgtraj_rast
    ForestChange <- subs(chgtraj_rast, forest_subs)
    if (cellStats(is.na(ForestChange), "sum") != ncell(ForestChange)) {
        plot_trajs(ForestChange, zoi,
                   data.frame(code=c(1, 2),
                              label=c("Forest gain", "Forest loss"),
                              color=c("#33FF33", "#FF3333"),
                              stringsAsFactors=FALSE),
                   title="Forest change")
        ggsave(paste0(file_path_sans_ext(chgtraj_file), '_transitions_natforest_change.png'), 
               width=img_width, height=img_height, dpi=img_dpi, type=type)
    }

    # Make plantation forest/non plantation forest transition image
    traj_codes$PlantChange <- NA
    traj_codes$PlantChange[with(traj_codes, (t0_name != "Plantation.forest") & (t1_name == "Plantation.forest"))] <- 1
    traj_codes$PlantChange[with(traj_codes, (t0_name == "Plantation.forest") & (t1_name != "Plantation.forest"))] <- 2
    planforest_subs <- data.frame(from=traj_codes$Code[!is.na(traj_codes$PlantChange)],
                                  to=traj_codes$PlantChange[!is.na(traj_codes$PlantChange)])
    PlantChange <- chgtraj_rast
    PlantChange <- subs(chgtraj_rast, planforest_subs)
    # Only make the plot if there are cells that show change
    if (cellStats(is.na(PlantChange), "sum") != ncell(PlantChange)) {
        plot_trajs(PlantChange, zoi,
                   data.frame(code=c(1, 2),
                              label=c("To plantation", "From plantation"),
                              color=c("#33FF33", "#FF3333"),
                              stringsAsFactors=FALSE),
                   title="Plantation change")
        ggsave(paste0(file_path_sans_ext(chgtraj_file), '_transitions_plantation_change.png'), 
               width=img_width, height=img_height, dpi=img_dpi, type=type)
    }

    # Make ag/non-ag transition image
    ag_classes <-  c("Agriculture", "Plantation.forest")
    traj_codes$AgChange <- NA
    traj_codes$AgChange[with(traj_codes, !(t0_name %in% ag_classes) & (t1_name %in% ag_classes))] <- 1
    traj_codes$AgChange[with(traj_codes, (t0_name %in% ag_classes) & !(t1_name %in% ag_classes))] <- 2
    ag_subs <- data.frame(from=traj_codes$Code[!is.na(traj_codes$AgChange)],
                              to=traj_codes$AgChange[!is.na(traj_codes$AgChange)])
    AgChange <- chgtraj_rast
    AgChange <- subs(chgtraj_rast, ag_subs)
    if (cellStats(is.na(AgChange), "sum") != ncell(AgChange)) {
        plot_trajs(AgChange, zoi,
                   data.frame(code=c(1, 2),
                              label=c("Ag gain", "Ag loss"),
                              color=c("#33FF33", "#FF3333"),
                              stringsAsFactors=FALSE),
                   title="Agricultural change")
        ggsave(paste0(file_path_sans_ext(chgtraj_file), '_transitions_ag_change.png'), 
               width=img_width, height=img_height, dpi=img_dpi, type=type)
    }

}

stopCluster(cl)
