# script to isolate example article

# for this reprex, we will use article_polygon A100_P1,
# which is a study of bird range shifts in Finland. 
# this script will isolate examples from that study to 
# run several processing steps that were completed to 
# derive values in the BioShifts database. 



# install BioShifts -------------------------------------------------------
# remotes::install_github("BioShifts/BioShiftR")
library(BioShiftR)
library(dplyr)
library(terra)
library(sf)
library(parallel)
library(pbapply)
library(rgbif)
library(data.table)
library(stringr)


# download study polygons (only needs to be done once):
download_polygons()


# isolate study -----------------------------------------------------------
# specify target article
target_article <- "A100"

# get all shifts
all_shifts <- get_shifts()

# filter to target
target_shifts <- all_shifts %>%
    filter(article_id == target_article)



## 1. estimate standardization  --------------------------------------------
# first, supplement target shifts with author reported rates and methods
# as provided in BioShiftR, which we will use to recalculate calc_rate.
standardization_demo <- target_shifts %>%
    add_author_reported() %>%
    add_methods()

source("workflows/01-estimate-standardization/01.1-standardization-function.R")
standardization_demo <- standardize_shifts(standardization_demo)

# check if demo calculations are same as in BioShiftR
table(standardization_demo$direction == standardization_demo$direction_demo)
table(round(standardization_demo$calc_rate,3) == round(standardization_demo$calc_rate_demo,3))
table(standardization_demo$calc_unit == standardization_demo$calc_unit_demo)

# remove demo df and function
rm(standardization_demo, standardize_shifts)


# 2. taxonomic harmonization ----------------------------------------------
taxo_demo <- target_shifts %>%
    distinct(sp_name_publication, sp_name_checked)

source("workflows/02-taxonomic-harmonization/02.1-taxonomic-harmonization-function.R")
taxo_demo <- harmonize_taxonomy(taxo_demo)

# check if demo taxonomies are the same as in BioShiftR
# NOTE that it is likely that some will differ as taxonomies continue to change
table(taxo_demo$sp_name_checked == taxo_demo$sp_name_checked_demo)

# remove demo and function
rm(taxo_demo, harmonize_taxonomy)


# 3. polygon summarization ------------------------------------------------

# download global elevation layer:
elevation <- geodata::elevation_global(res = "0.5", path = here::here("data-raw"))

# set equal area projection:
Moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
sf_use_s2(FALSE) 

# project to Mollweide (less deformation)
elevationM <- terra::project(
    elevation, Moll, threads = TRUE)


# isolate target polygon
demo_poly <- target_shifts %>%
    distinct(article_id, poly_id, type) %>%
    slice(1) %>%
    add_polygons() %>%
    add_poly_info()


# source the function
source("workflows/03-polygon-summarization/03.1-summarize-polygon-function.R")
demo_poly <- poly_sum(demo_poly)

demo_poly %>% glimpse()

# check if demo rates are same as in BioShiftR
setequal(demo_poly$lat_cent_deg, demo_poly$lat_cent_deg_demo)
setequal(demo_poly$lon_cent_deg, demo_poly$lon_cent_deg_demo)
setequal(demo_poly$lat_min, demo_poly$lat_min_demo)
setequal(demo_poly$lat_max, demo_poly$lat_max_demo)
setequal(demo_poly$lat_extent_km, demo_poly$lat_extent_km_demo)
setequal(demo_poly$area_km2, demo_poly$area_km2_demo)

# remove objects and function
rm(demo_poly, elevation, elevationM, Moll, poly_sum)



# 4. climate velocity -----------------------------------------------------

# isolate polygon to run CV (differs from above because includes duration)
# note that this polygon has 5 different durations. We will isolate 2 to 
# demonstrate a warming period and a cooling period.
cv_demo <- target_shifts %>%
    add_methods() %>%
    distinct(article_id, poly_id,eco, midpoint_firstperiod, midpoint_lastperiod, method_id, type) %>%
    add_cv(res = "25km", stat = c("q25","q75","median","mean","sd")) %>%
    add_trends(res = "25km", stat = c("mean","sd")) %>%
    add_baselines(res = "25km", stat = c("mean","sd")) %>%
    add_polygons() %>%
    filter( midpoint_firstperiod %in% c(1974.5, 1984.5),
            midpoint_lastperiod %in% c(1984.5, 1994.5)) 
# our reprex 

# settings 
res <- "25km" # only 25km data provided in this reprex

source("workflows/04-climate-velocity/04.1-cv-function.R")

# run demo 1 =======
cv_demo_1 <- calc_cv(cv_demo %>% slice(1))

# round values to 4 decimals
cv_demo_1 <- cv_demo_1 %>%
    mutate(across(where(is.numeric), \(x) round(x, digits = 4)))

# check cvs are the same:
cv_demo_1$cv_temp_mean - cv_demo_1$cv_lat.wtd.mean_demo
cv_demo_1$cv_temp_q25 - cv_demo_1$cv_lat.wtd.q25_demo
cv_demo_1$cv_temp_q75 - cv_demo_1$cv_lat.wtd.q75_demo
cv_demo_1$cv_temp_sd - cv_demo_1$cv_lat.wtd.stdev_demo

# check baselines are the same:
cv_demo_1$baseline_temp_mean - cv_demo_1$baseline.wtd.mean_demo
cv_demo_1$baseline_temp_sd - cv_demo_1$baseline.wtd.stdev_demo

# check trends are the same:
cv_demo_1$trend_temp_mean - cv_demo_1$ttrend.wtd.mean_demo
cv_demo_1$trend_temp_sd - cv_demo_1$ttrend.wtd.stdev_demo



# run demo 2 =======
cv_demo_2 <- calc_cv(cv_demo %>% slice(2))

# round values to 4 decimals
cv_demo_2 <- cv_demo_2 %>%
    mutate(across(where(is.numeric), \(x) round(x, digits = 4)))

# check cvs are the same:
cv_demo_2$cv_temp_mean - cv_demo_2$cv_lat.wtd.mean_demo
cv_demo_2$cv_temp_q25 - cv_demo_2$cv_lat.wtd.q25_demo
cv_demo_2$cv_temp_q75 - cv_demo_2$cv_lat.wtd.q75_demo
cv_demo_2$cv_temp_sd - cv_demo_2$cv_lat.wtd.stdev_demo

# check baselines are the same:
cv_demo_2$baseline_temp_mean - cv_demo_2$baseline.wtd.mean_demo
cv_demo_2$baseline_temp_sd - cv_demo_2$baseline.wtd.stdev_demo

# check trends are the same:
cv_demo_2$trend_temp_mean - cv_demo_2$ttrend.wtd.mean_demo
cv_demo_2$trend_temp_sd - cv_demo_2$ttrend.wtd.stdev_demo

# proportions
cv_demo_1$cv_temp_mean / cv_demo_1$cv_lat.wtd.mean_demo * 100
cv_demo_2$cv_temp_mean / cv_demo_2$cv_lat.wtd.mean_demo * 100

rm(list = ls())
