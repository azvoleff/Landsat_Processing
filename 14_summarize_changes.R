source('0_settings.R')

library(foreach)
library(itertools)

library(ggplot2)
library(gridExtra) # for unit
library(dplyr)
library(reshape2)
library(stringr)

img_height <- 10
img_width <- 10
img_dpi <- 300

overwrite <- TRUE

sites <- read.csv('Site_Code_Key.csv')
sitecodes <- sites$Site.Name.Code

out_dir <- 'C:/Users/azvoleff/Code/TEAM/teamlucc_scripts/Landsat_Processing/LUCC_Results'
stopifnot(file_test('-d', out_dir))

for (sitecode in sitecodes) {
    message(paste0('Calculating stats for ', sitecode, '...'))
    image_basedir <- file.path(prefix, 'Landsat', sitecode)
    freqs_files <- dir(image_basedir,
                       pattern='^[a-zA-Z]*_[0-9]{4}-[0-9]{4}_chgdetect_chgtraj_freqs.csv$')

    if (length(freqs_files) == 0) {
        next
    }

    class_key <- read.csv('H:/Data/TEAM/LCLUC_Training/LCLUC_Classes_Key.csv')

    sitename <- as.character(sites$Site.Name.Short[match(sitecode, sites$Site.Name.Code)])

    foreach(freqs_file=iter(freqs_files),
            .packages=c('ggplot2', 'dplyr')) %do% {

        time_string <- str_extract(freqs_file, '[0-9]{4}-[0-9]{4}')

        freqs <- read.csv(file.path(image_basedir, freqs_file))
        freqs$t0_name <- gsub('X', '', freqs$t0_name)
        freqs$t0_name <- class_key$Class[match(freqs$t0_name, class_key$Code)]
        freqs$t1_name <- gsub('X', '', freqs$t1_name)
        freqs$t1_name <- class_key$Class[match(freqs$t1_name, class_key$Code)]
        freqs$t0_name_abbrev <- abbreviate(freqs$t0_name, 6)
        freqs$t1_name_abbrev <- abbreviate(freqs$t1_name, 6)
        freqs$Transition <- paste(freqs$t0_name_abbrev, freqs$t1_name_abbrev, sep=' -> ')

        classes2plot <- c('Plantation forest', 'Natural forest', 'Bare')

        for (class2plot in classes2plot) {
            ggplot(freqs[freqs$t0_name == class2plot, ]) +
                geom_bar(aes(Transition, freq), stat='identity') +
                ggtitle(paste(sitename, '-', time_string))
            ggsave(file.path(out_dir, paste0('transitions_', sitecode, 
                                             '_barplot_', class2plot, '_', 
                                             time_string, '.png')),
                   height=img_height, width=img_width, dpi=img_dpi)
        }


        # trans_mat <- dcast(data.frame(t0=freqs$t0_name,
        #                               t1=freqs$t1_name,
        #                               freq=freqs$freq),
        #                    t0 ~ t1, value.var='freq')

        freqs$t0_name <- ordered(freqs$t0_name, levels=rev(c('Urban/built',
                                                             'Agriculture',
                                                             'Plantation forest',
                                                             'Natural forest',
                                                             'Other vegetation',
                                                             'Bare',
                                                             'Water')))
        freqs$t1_name <- ordered(freqs$t1_name,
                                 levels=rev(levels(freqs$t0_name)))

        # TODO: Need to add in classes with zero occurrence to ensure all 
        # factors are represented in the freqs list (to ensure that all the 
        # plots have the same rows/columns across sites

        persist <- summarize(group_by(freqs, t0_name),
                             pct=paste0(round(freq[t0_name == t1_name] / sum(freq), 2)))

        percent_persistence <- freqs[freqs$t0_name == freqs$t1_name, ]

        freqs_zeroed_persist <- freqs
        freqs_zeroed_persist$freq[freqs_zeroed_persist$t0_name == freqs_zeroed_persist$t1_name] <- 0
        ggplot(freqs_zeroed_persist) +
            theme_bw() +
            geom_tile(aes(x=t1_name, y=t0_name, fill=freq/sum(freq)), colour='black') +
            geom_text(aes(x=t0_name, y=t0_name, label=pct), data=persist) +
            scale_fill_gradientn('Relative\nFrequency', limits=c(0, .25),
                                 colours=c('white', 'orange', 'red')) +
            xlab(paste('Class in', str_extract(time_string, '[0-9]{4}$'))) +
            ylab(paste('Class in', str_extract(time_string, '^[0-9]{4}'))) +
            theme(axis.text.y=element_text(angle=90, hjust=.5),
                  legend.key.size=unit(1.5, "line"),
                  panel.grid.major=element_blank()) +
        ggsave(file.path(out_dir, paste0('transitions_', sitecode, 
                                               '_colorplot_', time_string, 
                                               '.png')),
               height=img_height, width=img_width, dpi=img_dpi)
    }
}
