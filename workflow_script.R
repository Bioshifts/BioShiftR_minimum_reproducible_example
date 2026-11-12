# script to isolate example article

# for this reprex, we will use article_polygon A100_P1,
# which 

remove.packages("BioShiftR")
remotes::install_github("BioShifts/BioShiftR", force = T)

library(BioShiftR)
library(dplyr)

# download polygons locally (only needs to be done once)
# download_polygons()

# get all shifts
all_shifts <- get_shifts()

# specift target article
target_article <- "A100"

# isolate polygon to run polygon summarization
target_poly <- all_shifts %>%
    filter(article_id == target_article) %>%
    distinct(article_id, poly_id) %>%
    add_polygons() 

# isolate polygon to run CV (differs from above because includes duration)
# note that this polygon has 5 different durations. We will isolate 2 to 
# demonstrate a warming period and a cooling period.
target_polygon_cv <- all_shifts %>%
    filter(article_id == target_article) %>%
    add_methods() %>%
    distinct(article_id, poly_id, midpoint_firstperiod, midpoint_lastperiod) %>%
    add_polygons() %>%
    filter( midpoint_firstperiod %in% c(1974.5, 1984.5),
            midpoint_lastperiod %in% c(1984.5, 1994.5)) 


# isolate species names from publication for harmonization
target_spp_publication <- all_shifts %>%
    filter(article_id == target_article) %>%
    distinct(sp_name_publication) %>%
    # optional slice to isolate one:
    #slice(1) %>%
    pull()

# view how many names should change in harmonization
all_shifts %>%
    filter(article_id == target_article) %>%
    distinct(sp_name_checked, sp_name_publication) %>%
    mutate(same_name = sp_name_checked == sp_name_publication) %>%
    count(same_name)
# so 20 out of 94 should change in harmonization

# isolate author reported rates for standardization
# (in this case, they just report a poleward distance moved, so it's easy)
target_author_reported <- 
    all_shifts %>%
    filter(article_id == target_article) %>%
    add_author_reported() %>% 
    add_methods() %>% 
    select(sp_name_publication, 
           calc_rate,
           contains("author_reported"),
           contains("period"),
           duration) 

# count how many the freshly calculated rate matches calc_rate
target_author_reported %>% 
    mutate(test = as.numeric(author_reported)/duration) %>% 
    count(round(calc_rate,2) == round(test,2))


