source('0_settings.R')

espa_email <- 'azvoleff@gmail.com'

###############################################################################
# All Sites
download_folder <- 'I:/Landsat_Originals'

stopifnot(file_test('-d', download_folder))

# espa_download(espa_email, '5142014-155630', download_folder) # 0-2% cover
# espa_download(espa_email, '5142014-155558', download_folder) # 2-5% cover
# espa_download(espa_email, '5142014-15558', download_folder) # 5-10% cover
# espa_download(espa_email, '5142014-153954', download_folder) # 10-20% cover
espa_download(espa_email, '6162014-8585', download_folder) # 10-20% cover
