# Install packages #####

#install.packages("devtools")
library(tidyverse)
library(devtools)
library(here)
devtools::install_github("https://github.com/MacroSHEDS/macrosheds.git")
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

# Select target solute #####
target_solute = 'SO4_S'
# Read in data from the larger dataset #####
target_domain = 'hbef'
q <- ms_load_product(
    macrosheds_root = macrosheds_root,
    prodname = 'discharge',
    #filter_vars,
    #networks = 'lter',
    domains = target_domain,
    #site_codes = 'all',
    sort_result = FALSE
    )

chem <- ms_load_product(
    macrosheds_root = macrosheds_root,
    prodname = 'stream_chemistry',
    filter_vars = 'SO4_S',
    #networks = 'lter',
    domains = target_domain,
    warn = F
)

# select one site
target_site <- 'w1'

flux <- read_feather(here('ms_flux_annual', target_domain,
                          'stream_flux', paste0(target_site, '.feather'))) %>%
    filter(var == paste0('GN_', target_solute),
           ms_recommended == 1)

model <- summary(lm(data = flux, val~as.integer(wy), na.action = NULL))


ggplot(flux, aes(x = as.integer(wy), y = val))+
    geom_point()+
    geom_smooth(method = 'lm')


