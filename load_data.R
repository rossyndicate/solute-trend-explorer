# Intro #####
# Purpose is to load data in for solute trend detection. Data is in MacroSheds and
# on figshare
# Install packages #####
#install.packages("devtools")
library(tidyverse)
library(devtools)
library(here)
#devtools::install_github("https://github.com/MacroSHEDS/macrosheds.git")
library(macrosheds)
library(rfigshare)
library(feather)

# Set up MS package and load data from it #####
macrosheds_root <- here('ms_data')
ms_download_core_data(macrosheds_root = here('ms_data'),
                      networks = 'lter')

ms_download_ws_attr(macrosheds_root = macrosheds_root)
summaries <- read_feather(here('ms_data', 'watershed_summaries.feather'))

# uncomment below to load ms data (takes a while)
#ms_download_core_data(macrosheds_root = macrosheds_root, networks = 'lter', quiet = FALSE)

macrosheds::ms_load_sites()

# Read in flux data from figshare #####
info <- fs_download(article_id = 24975504, mine = F, session = NULL)
download.file(info, destfile = 'annual_loads.zip', mode="wb")
unzip('annual_loads.zip')




