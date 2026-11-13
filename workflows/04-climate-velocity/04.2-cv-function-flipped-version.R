# Script to Demo Latitudinal Climate Velocity Calculation



# libraries ---------------------------------------------------------------
library(BioShiftR)
library(terra)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyterra)
library(sf)
library(fields)
library(exactextractr)
theme_set(theme_bw())


# source climate velocity functions ---------------------------------------
source("workflows/04-climate-velocity/04.0-cv-internals.R")


calc_cv_2 <- function(data){
    
    
    if(nrow(data) > 1)stop("this reprex runs one polygon at a time")
    
    # change start/end dates to match data availability -------------------------
    data <- data %>%
        mutate(
            midpoint_firstperiod_cv = case_when(
                eco == "Ter" & midpoint_firstperiod >= 1902 ~ midpoint_firstperiod,
                eco == "Ter" & midpoint_firstperiod < 1902 ~ 1902,
                eco == "Mar" & midpoint_firstperiod >= 1959 ~ midpoint_firstperiod,
                eco == "Mar" & midpoint_firstperiod < 1959 ~ 1959
            )
        )
    
    # isolate polygon
    poly <- vect(data)
    
    # get focal years
    start_year <- floor(poly$midpoint_firstperiod_cv) # round down for start year
    end_year <- ceiling(poly$midpoint_lastperiod) # round up for end year
    eco <- poly$eco # Mar or Ter; designates which temperature dataset to use
    
    # list temperature years to upload
    upload_years <- start_year:end_year
    # make vector of file names to upload
    if(eco == "Ter") {
        upload_files <- paste0("bios_Ter_",upload_years,".tif")
    }
    if(eco == "Mar"){
        upload_files <-paste0("bios_Mar_",upload_years,".tif")
    }
    
    # upload rasters 
    temp_data_path <- "data-raw/temperature-data"
    path <- file.path(temp_data_path, 
                      paste0("res-",res),
                      ifelse(eco == "Ter","terrestrial","marine"))
    
    # upload and rasterize files
    r <- rast(file.path(path,upload_files))
    # add year to layer names
    names(r) <- paste0("mat_", upload_years)
    r <- r * 1
    # plot to see a layer:
    # plot(r[[1]])
    
    # 1. crop to buffered polygon ---------------------------------------------
    # create polygon buffer around the study polygon to incorporate cells from approximately
    # 1 cell width outside the polygon into calculations 
    buffer_size_m <- as.numeric(str_remove(res,"km"))*1000 
    buffered_poly <- terra::buffer(poly, width = buffer_size_m)
    
    # crop raster
    terra::window(r) <- terra::ext(buffered_poly)
    r <- terra::mask(r, buffered_poly)
    terra::window(r) <- NULL
    
    # convert to real degrees (chelsa uses degrees*10)
    if(eco == "Ter") r <- r/10 
    gc() # clear memory after reducing raster size
    # plot to view:
    #ggplot() +
    #    geom_spatraster(data = r[[c(1,round(nlyr(r)/2),nlyr(r))]]) +
    #    facet_wrap(~lyr) +
    #    scale_fill_viridis_c(na.value = "transparent") +
    #    geom_spatvector(data = poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1) + 
    #    geom_spatvector(data = buffered_poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1,
    #                    linetype = "dashed")
    
    
    # 2.1 find temperature trend ----------------------------------------------
    # `temp_grad` applies a custom function (`temp_gradFun`) to every raster cell,
    # resulting in a raster layer of the warming trend of each cell (the slope of the linear 
    # model of temperature in each cell over years).
    ttrend = temp_trend(
        r,
        th = 0.25*nlyr(r)) ## set minimum N obs. to 1/4 time series length
    # plot to view:
    #ggplot() +
    #    geom_spatraster(data = ttrend) +
    #    scale_fill_gradient2(low = "cornflowerblue", high = "tomato2", na.value = "transparent") +
    #    geom_spatvector(data = poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1) + 
    #    geom_spatvector(data = buffered_poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1,
    #                    linetype = "dashed")
    
    # 2.2 find spat grad -------------------------------------------------
    # find average temperature in each cell:
    avg_r <- terra::app(r, fun = mean, na.rm = TRUE)
    
    # plot to view:
    #ggplot() +
    #    geom_spatraster(data = avg_r) +
    #    scale_fill_viridis_c(na.value = "transparent") +
    #    geom_spatvector(data = poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1) + 
    #    geom_spatvector(data = buffered_poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1,
    #                    linetype = "dashed")
    
    # calculate spatial gradient of mean temperature 
    # (the angle and magnitude of average warming)
    spgrad = spatial_grad_JL(avg_r)
    
    # flip the spatial gradient so it's pointing towards cold instead of warm
    spgrad_flipped <- spgrad_flip(spgrad)
    rm(spgrad)
    
    # 3.1 find directional climate velocity -----------------------------------
    if(ncell(spgrad_flipped) > 100){ truncate = T} else {truncate = F}
    CV_flipped <- gVelocity_JL_flipped(grad = spgrad_flipped, slope = ttrend, truncate = truncate)
    
    # 3.2 Find latitudinal CV --------------------------------------------------
    # extract the northward vector component:
    CV_lat <- CV_flipped$Vel * cos(deg_to_rad(CV_flipped$Ang))
    # crop back to cells within the focal polygon (removing the buffer cells):
    CV_lat <- mask(CV_lat, poly)
    
    # 3.3 Flip sign for southern hemispherre -----------------------------------
    # if polygon is in southern hemisphere, flip the sign of lat velocity
    # so that it becomes "poleward" instead of "northward" velocity
    lat <- terra::init(CV_lat, "y") |> terra::mask(CV_lat)
    min_lat <- as.numeric(global(lat,min,na.rm=TRUE))
    if(min_lat < 0){
        CV_lat <- ifel(lat < 0, CV_lat * -1, CV_lat)
    }
    rm(lat, min_lat)
    
    # plot CV and Latitudinal CV
    #cv_checks_2 <- jl_plot3_gvelLat(rast_gvel = CV, 
    #                 rast_gvelLat = CV_lat, 
    #                 avg_rast = avg_r,
    #                 SA = poly,
    #                 spoke_fact = 4,
    #                 gradfact = 1
    #                 )
    #cv_checks_2
    # here, black arrows represent the "directional" cv (at some angle), 
    # and red arrows represent the vector component of the angle facing 0°N
    
    # 4. find average inside polygon ------------------------------------------
    # calculate a few metrics of CV and variance, weighted by cell area and coverage fraction
    
    # get all values of raster cells
    vals_cv <- values(CV_lat)[,1] #  latitudinal CVs
    vals_ttrend <- values(ttrend)[,1] # temperature trends
    vals_baseline <- values(avg_r)[,1]
    # make weights vector: coverage fraction of cells * cell size
    weights <- values(
        exactextractr::coverage_fraction(CV_lat, st_as_sf(poly))[[1]] * 
            cellSize(CV_lat, mask = T)
    )[,1]
    
    # calculate all weighted values:
    wt.mean <- Hmisc::wtd.mean(vals_cv, weights)
    wt.stdev <- sqrt(Hmisc::wtd.var(vals_cv, weights))
    wt.quants <- Hmisc::wtd.quantile(vals_cv, weights, probs = c(.25,.5,.75)) %>% unname()
    wt.mean.ttrend <- Hmisc::wtd.mean(vals_ttrend, weights)
    wt.stdev.ttrend <- sqrt(Hmisc::wtd.var(vals_ttrend, weights))
    wt.mean.baseline <- Hmisc::wtd.mean(vals_baseline, weights)
    wt.stdev.baseline <- sqrt(Hmisc::wtd.var(vals_baseline, weights))
    
    out.vals <- data.frame(
        "cv_lat.wtd.mean" = wt.mean,
        "cv_lat.wtd.stdev" = wt.stdev,
        "cv_lat.wtd.q25" = wt.quants[1],
        "cv_lat.wtd.q50" = wt.quants[2],
        "cv_lat.wtd.q75" = wt.quants[3],
        "ttrend.wtd.mean" = wt.mean.ttrend,
        "ttrend.wtd.stdev" = wt.stdev.ttrend,
        "baseline.wtd.mean" = wt.mean.baseline,
        "baseline.wtd.stdev" = wt.stdev.baseline
    )
    rm(wt.mean, wt.stdev, wt.quants, wt.mean.ttrend, wt.stdev.ttrend, wt.mean.baseline, wt.stdev.baseline)
    
    colnames(out.vals) <- paste0(colnames(out.vals),"_demo")
    out.vals$article_id <- data$article_id
    out.vals$poly_id <- data$poly_id
    
    data <- data %>%
        left_join(out.vals)
    
    
}


