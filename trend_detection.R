library(tidyverse)
library(feather)
library(macrosheds)
library(lubridate)
macrosheds_root <- here('ms_data')


# Select targets #####
target_solute = 'SO4_S'
target_domain = 'hbef'
site_list <- tools::file_path_sans_ext(list.files(here('ms_flux_annual',
                                                       target_domain, 'stream_flux'),
                                       pattern = '.feather$'))

out_frame <- tibble(domain = as.character(),
                    site = as.character(),
                    solute = as.character(),
                    n = as.integer(),
                    con_trend_slope = as.numeric(),
                    con_trend_p = as.numeric(),
                    con_trend_r_squared = as.numeric(),
                    flux_trend_slope = as.numeric(),
                    flux_trend_p = as.numeric(),
                    flux_trend_r_squared = as.numeric())

for(i in 1:length(site_list)){
target_site <- site_list[i]
# load in flux data #####

flux <- read_feather(here('ms_flux_annual', target_domain,
                          'stream_flux', paste0(target_site, '.feather'))) %>%
    filter(var == paste0('GN_', target_solute),
           ms_recommended == 1)
    if(nrow(flux>14)){
# Read in data from the larger dataset #####
q <- ms_load_product(
    macrosheds_root = macrosheds_root,
    prodname = 'discharge',
    #filter_vars,
    #networks = 'lter',
    domains = target_domain,
    site_codes = site_list,
    sort_result = FALSE
)

chem <- ms_load_product(
    macrosheds_root = macrosheds_root,
    prodname = 'stream_chemistry',
    filter_vars = target_solute,
    #networks = 'lter',
    domains = target_domain,
    warn = F,
    site_codes = site_list
) %>%
    mutate(wy = as.integer(year(datetime))) %>%
    filter(!is.na(val)) %>%
    group_by(wy) %>%
    summarize(mean_con = mean(val))

con_model <- summary(lm(data = chem, mean_con~wy, na.action = NULL))

# make flux model #####
flux_model <- summary(lm(data = flux, val~as.integer(wy), na.action = NULL))

# add outputs to table #####
loop_frame <- tibble(domain = target_domain,
       site = target_site,
       solute = target_solute,
       n = nrow(flux),
       con_trend_slope = con_model$coefficients[2,1],
       con_trend_p = con_model$coefficients[2,4],
       con_trend_r_squared = con_model$r.squared,
       flux_trend_slope = flux_model$coefficients[2,1],
       flux_trend_p = flux_model$coefficients[2,4] ,
       flux_trend_r_squared = flux_model$r.squared)
out_frame <- rbind(out_frame, loop_frame)
    } #end data checking if/else
else{next}
}# end loop
