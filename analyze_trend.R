library(tidyverse)
library(feather)
library(macrosheds)
library(lubridate)
library(here)
library(lfstat)
macrosheds_root <- here('ms_data')

#functions for later #####
aggregate_seasonal_data <- function(target_season){
    con_sea <- con_tag %>%
        dplyr::filter(season == target_season) %>%
        group_by(season_year) %>%
        summarize(mean_con = mean(val),
                  high_con = as.numeric(quantile(val, 0.95)),
                  low_con = as.numeric(quantile(val, 0.05)),
                  date_mean = as.Date(mean(datetime)))
}

# initialize outputs ####
out_frame <- tibble(
    # top level
    domain = as.character(),
    site = as.character(),
    solute = as.character(),

    # annual
    n_ann = as.integer(),
    slope_ann_mean = as.numeric(),
    p_ann_mean = as.numeric(),
    rsquared_ann_mean = as.numeric(),

    slope_ann_low = as.numeric(),
    p_ann_low = as.numeric(),
    rsquared_ann_low = as.numeric(),

    slope_ann_high = as.numeric(),
    p_ann_high = as.numeric(),
    rsquared_ann_high = as.numeric(),

    # spring
    slope_sp_mean = as.numeric(),
    p_sp_mean = as.numeric(),
    rsquared_sp_mean = as.numeric(),

    slope_sp_low = as.numeric(),
    p_sp_low = as.numeric(),
    rsquared_sp_low = as.numeric(),

    slope_sp_high = as.numeric(),
    p_sp_high = as.numeric(),
    rsquared_sp_high = as.numeric(),

    # winter
    slope_w_mean = as.numeric(),
    p_w_mean = as.numeric(),
    rsquared_w_mean = as.numeric(),

    slope_w_low = as.numeric(),
    p_w_low = as.numeric(),
    rsquared_w_low = as.numeric(),

    slope_w_high = as.numeric(),
    p_w_high = as.numeric(),
    rsquared_w_high = as.numeric(),

    # summer
    slope_su_mean =as.numeric(),
    p_su_mean = as.numeric(),
    rsquared_su_mean = as.numeric(),

    slope_su_low = as.numeric(),
    p_su_low = as.numeric(),
    rsquared_su_low = as.numeric(),

    slope_su_high = as.numeric(),
    p_su_high = as.numeric(),
    rsquared_su_high = as.numeric(),

    # fall
    slope_f_mean = as.numeric(),
    p_f_mean = as.numeric(),
    rsquared_f_mean = as.numeric(),

    slope_f_low = as.numeric(),
    p_f_low = as.numeric(),
    rsquared_f_low = as.numeric(),

    slope_f_high = as.numeric(),
    p_f_high = as.numeric(),
    rsquared_f_high = as.numeric())


# Select targets #####
target_solute = 'SO4_S'

# find sites with data ####

site_data <- ms_load_sites()
domain_list <- list.dirs(here('ms_data'), recursive = F, full.names = F)

#remove non-conus domains if present
domain_list <- domain_list[!domain_list %in% c('arctic', 'mcmurdo', 'krycklan', 'luquillo')]

# loop through domains
for(u in 1:length(domain_list)){
target_domain = domain_list[u]

site_list <- tools::file_path_sans_ext(list.files(here('ms_data',
                                                       target_domain, 'stream_chemistry'),
                                       pattern = '.feather$'))
if(length(site_list > 0)){

# loop through sites ####
    for(i in 1:length(site_list)){
    target_site <- site_list[i]
    # load in flux data #####

    con <- read_feather(here('ms_data', target_domain,
                              'stream_chemistry', paste0(target_site, '.feather'))) %>%
        filter(var == paste0('GN_', target_solute)) %>%
        select(datetime, val) %>%
        na.omit()
        if(nrow(con > 120)){
    # sites past this point will have enough data to make ten years with seasons #####
    # now make sure there is good data coverage ####
    con_tag <- con %>%
        mutate(wy = water_year(datetime, origin = 'usgs'),
               year = year(datetime),
               month = month(datetime),
               month_year = paste0(month,'_',year),
               season = '')

    # set seasons
    con_tag$season[con_tag$month %in% c(12,1,2)] <- 'w'
    con_tag$season[con_tag$month %in% c(3,4,5)] <- 'sp'
    con_tag$season[con_tag$month %in% c(6,7,8)] <- 'su'
    con_tag$season[con_tag$month %in% c(9,10,11)] <- 'f'

    con_tag <- con_tag %>%
        mutate(season_year = paste0(season,'_',year))

    # check for at least 12 samples in each year
    year_check <- con_tag %>%
        group_by(wy) %>%
        summarize(n = n()) %>%
        filter(n > 11)

    if(nrow(year_check) > 9){test_years = 1}else{test_years = 0}

    # do the same thing for seasons, need at least 4 samples in a season to make a fit
    # so need min of 40 seasons for 10 years of data
    sea_check <- con_tag %>%
        group_by(season_year) %>%
        summarize(n = n(),
                  year = mean(year),
                  season = season) %>%
        filter(n > 3)

    if(nrow(sea_check) > 9){test_seas_long = 1}else{test_seas_long = 0}
    # also want to make sure they sample across at least 3 of the four seasons
    if(length(unique(sea_check$season)) > 2){test_seas_even = 1}else{test_seas_even = 0}

    # if any checks fail, skip site and move to the next #####
    if(test_seas_long*test_seas_even*test_years==1){
    # Read in data #####
    chem_ann <- con_tag %>%
        group_by(wy = as.integer(as.character(wy))) %>%
        summarize(mean_con = mean(val),
                  high_con = as.numeric(quantile(val, 0.95)),
                  low_con = as.numeric(quantile(val, 0.05)))

    # make models of annual data #####
    ann_model_mean <- summary(lm(data = chem_ann, mean_con~wy, na.action = NULL))
    ann_model_high <- summary(lm(data = chem_ann, high_con~wy, na.action = NULL))
    ann_model_low <- summary(lm(data = chem_ann, low_con~wy, na.action = NULL))

    # do the same for each season #####
    spring_data <- aggregate_seasonal_data('sp')
    sp_model_mean <- summary(lm(data = spring_data, mean_con~date_mean, na.action = NULL))
    sp_model_high <- summary(lm(data = spring_data, high_con~date_mean, na.action = NULL))
    sp_model_low <- summary(lm(data = spring_data, low_con~date_mean, na.action = NULL))

    summer_data <- aggregate_seasonal_data('su')
    su_model_mean <- summary(lm(data = summer_data, mean_con~date_mean, na.action = NULL))
    su_model_high <- summary(lm(data = summer_data, high_con~date_mean, na.action = NULL))
    su_model_low <- summary(lm(data = summer_data, low_con~date_mean, na.action = NULL))

    fall_data <- aggregate_seasonal_data('f')
    f_model_mean <- summary(lm(data = fall_data, mean_con~date_mean, na.action = NULL))
    f_model_high <- summary(lm(data = fall_data, high_con~date_mean, na.action = NULL))
    f_model_low <- summary(lm(data = fall_data, low_con~date_mean, na.action = NULL))


    winter_data <- aggregate_seasonal_data('w')
    if(nrow(winter_data) < 3){
        no_winter <- 1 #handle no winter data
        next
    }else{
    no_winter = 0 # handle no winter data
    w_model_mean <- summary(lm(data = winter_data, mean_con~date_mean, na.action = NULL))
    w_model_high <- summary(lm(data = winter_data, high_con~date_mean, na.action = NULL))
    w_model_low <- summary(lm(data = winter_data, low_con~date_mean, na.action = NULL))
    }



    # add outputs to table #####
    if(no_winter == 0){
    loop_frame <- tibble(
           # top level
           domain = target_domain,
           site = target_site,
           solute = target_solute,

           # annual
           n_ann = nrow(chem_ann),
           slope_ann_mean = ann_model_mean$coefficients[2,1],
           p_ann_mean = ann_model_mean$coefficients[2,4],
           rsquared_ann_mean = ann_model_mean$r.squared,

           slope_ann_low = ann_model_low$coefficients[2,1],
           p_ann_low = ann_model_low$coefficients[2,4],
           rsquared_ann_low = ann_model_low$r.squared,

           slope_ann_high = ann_model_high$coefficients[2,1],
           p_ann_high = ann_model_high$coefficients[2,4],
           rsquared_ann_high = ann_model_high$r.squared,

           # spring
           slope_sp_mean = sp_model_mean$coefficients[2,1],
           p_sp_mean = sp_model_mean$coefficients[2,4],
           rsquared_sp_mean = sp_model_mean$r.squared,

           slope_sp_low = sp_model_low$coefficients[2,1],
           p_sp_low = sp_model_low$coefficients[2,4],
           rsquared_sp_low = sp_model_low$r.squared,

           slope_sp_high = sp_model_high$coefficients[2,1],
           p_sp_high = sp_model_high$coefficients[2,4],
           rsquared_sp_high = sp_model_high$r.squared,

           # winter
           slope_w_mean = w_model_mean$coefficients[2,1],
           p_w_mean = w_model_mean$coefficients[2,4],
           rsquared_w_mean = w_model_mean$r.squared,

           slope_w_low = w_model_low$coefficients[2,1],
           p_w_low = w_model_low$coefficients[2,4],
           rsquared_w_low = w_model_low$r.squared,

           slope_w_high = w_model_high$coefficients[2,1],
           p_w_high = w_model_high$coefficients[2,4],
           rsquared_w_high = w_model_high$r.squared,

           # summer
           slope_su_mean = su_model_mean$coefficients[2,1],
           p_su_mean = su_model_mean$coefficients[2,4],
           rsquared_su_mean = su_model_mean$r.squared,

           slope_su_low = su_model_low$coefficients[2,1],
           p_su_low = su_model_low$coefficients[2,4],
           rsquared_su_low = su_model_low$r.squared,

           slope_su_high = su_model_high$coefficients[2,1],
           p_su_high = su_model_high$coefficients[2,4],
           rsquared_su_high = su_model_high$r.squared,

           # fall
           slope_f_mean = f_model_mean$coefficients[2,1],
           p_f_mean = f_model_mean$coefficients[2,4],
           rsquared_f_mean = f_model_mean$r.squared,

           slope_f_low = f_model_low$coefficients[2,1],
           p_f_low = f_model_low$coefficients[2,4],
           rsquared_f_low = f_model_low$r.squared,

           slope_f_high = f_model_high$coefficients[2,1],
           p_f_high = f_model_high$coefficients[2,4],
           rsquared_f_high = f_model_high$r.squared)} # end yes winter init
    else{
        loop_frame <- tibble(
            # top level
            domain = target_domain,
            site = target_site,
            solute = target_solute,

            # annual
            n_ann = nrow(chem_ann),
            slope_ann_mean = ann_model_mean$coefficients[2,1],
            p_ann_mean = ann_model_mean$coefficients[2,4],
            rsquared_ann_mean = ann_model_mean$r.squared,

            slope_ann_low = ann_model_low$coefficients[2,1],
            p_ann_low = ann_model_low$coefficients[2,4],
            rsquared_ann_low = ann_model_low$r.squared,

            slope_ann_high = ann_model_high$coefficients[2,1],
            p_ann_high = ann_model_high$coefficients[2,4],
            rsquared_ann_high = ann_model_high$r.squared,

            # spring
            slope_sp_mean = sp_model_mean$coefficients[2,1],
            p_sp_mean = sp_model_mean$coefficients[2,4],
            rsquared_sp_mean = sp_model_mean$r.squared,

            slope_sp_low = sp_model_low$coefficients[2,1],
            p_sp_low = sp_model_low$coefficients[2,4],
            rsquared_sp_low = sp_model_low$r.squared,

            slope_sp_high = sp_model_high$coefficients[2,1],
            p_sp_high = sp_model_high$coefficients[2,4],
            rsquared_sp_high = sp_model_high$r.squared,

            # winter
            slope_w_mean = NA,
            p_w_mean = NA,
            rsquared_w_mean = NA,

            slope_w_low = NA,
            p_w_low = NA,
            rsquared_w_low = NA,

            slope_w_high = NA,
            p_w_high = NA,
            rsquared_w_high = NA,

            # summer
            slope_su_mean = su_model_mean$coefficients[2,1],
            p_su_mean = su_model_mean$coefficients[2,4],
            rsquared_su_mean = su_model_mean$r.squared,

            slope_su_low = su_model_low$coefficients[2,1],
            p_su_low = su_model_low$coefficients[2,4],
            rsquared_su_low = su_model_low$r.squared,

            slope_su_high = su_model_high$coefficients[2,1],
            p_su_high = su_model_high$coefficients[2,4],
            rsquared_su_high = su_model_high$r.squared,

            # fall
            slope_f_mean = f_model_mean$coefficients[2,1],
            p_f_mean = f_model_mean$coefficients[2,4],
            rsquared_f_mean = f_model_mean$r.squared,

            slope_f_low = f_model_low$coefficients[2,1],
            p_f_low = f_model_low$coefficients[2,4],
            rsquared_f_low = f_model_low$r.squared,

            slope_f_high = f_model_high$coefficients[2,1],
            p_f_high = f_model_high$coefficients[2,4],
            rsquared_f_high = f_model_high$r.squared)}# end no winter output int

    out_frame <- rbind(out_frame, loop_frame)


    } # end if check on years and season
    else{next}
        } #end if check on site level data length
    else{next}
    }# end site loop
}# end site list length check
else{next}
}# end domain loop
