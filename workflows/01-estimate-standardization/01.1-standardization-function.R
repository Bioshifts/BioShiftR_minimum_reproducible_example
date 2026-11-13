

# 1. Standardization

# Shifts are reported from source articles in various formats (distance that a
# range moved over time, velocity in some unit, distance at an angle, etc). 
# This script demonstrates how author reported rates are standardized towards
# km/y in the poleward direction for latitudinal shifts 
# and m/year upwards for elevation shifts


standardize_shifts <- function(data){
    
    # convert author reported to numeric
    data$author_reported <- suppressWarnings(as.numeric(data$author_reported))
    
    # calculate the shift direction
    data <- data %>%
        mutate(
            direction_demo = case_when(
                author_reported < 0 & param == "TE" ~ "Expansion",
                author_reported > 0 & param == "TE" ~ "Contraction",
                author_reported < 0 & param == "LE" ~ "Contraction",
                author_reported > 0 & param == "LE" ~ "Expansion",
                author_reported > 0 & param == "O" & type == "ELE" ~ "Higher Elevation",
                author_reported > 0 & param == "O" & type == "LAT" ~ "Towards Poles",
                author_reported < 0 & param == "O" & type == "ELE" ~ "Lower Elevation",
                author_reported < 0 & param == "O" & type == "LAT" ~ "Towards Equator",
                author_reported %in% c("0", "0.0", "NULL") ~ "No Change",
                TRUE ~ NA_character_
            )
        )
        
    # calculate estimated shift rates
    KM_PER_DEG_LAT <- 111.32
    
    
    data <- data %>%
        mutate(
            Rate_v2 = case_when(
                author_reported_unit == "deg"       ~ author_reported * KM_PER_DEG_LAT / duration,
                author_reported_unit == "deg/year"  ~ author_reported * KM_PER_DEG_LAT,
                author_reported_unit == "km"        ~ author_reported / duration,
                author_reported_unit == "m"         ~ author_reported / duration,
                author_reported_unit == "m/decade"  ~ author_reported / 10,
                author_reported_unit == "km/decade" ~ author_reported / 10,
                TRUE ~ author_reported
            ),
            UnitRate_v2 = case_when(
                author_reported_unit %in% c("deg","deg/year","km") ~ "km/year",
                author_reported_unit %in% c("m","m/decade")        ~ "m/year",
                author_reported_unit == "km/decade"                ~ "km/year",
                TRUE ~ author_reported_unit
            )
        )
    
    
    # unify elevation records to m/year and latitude records to km/year
    data <- data %>%
        mutate(
            Rate_v2 = case_when(
                type == "ELE" & UnitRate_v2 == "km/year" ~ Rate_v2 * 1000,
                type == "LAT" & UnitRate_v2 == "m/year"  ~ Rate_v2 / 1000,
                TRUE ~ Rate_v2
            ),
            UnitRate_v2 = case_when(
                type == "ELE" & UnitRate_v2 == "km/year" ~ "m/year",
                type == "LAT" & UnitRate_v2 == "m/year"  ~ "km/year",
                TRUE ~ UnitRate_v2
            )
        )
    
    data$Rate_v2 <- round(data$Rate_v2,4)
    
    data <- data %>%
        rename(calc_rate_demo = Rate_v2,
               calc_unit_demo = UnitRate_v2)

    return(data)
}





