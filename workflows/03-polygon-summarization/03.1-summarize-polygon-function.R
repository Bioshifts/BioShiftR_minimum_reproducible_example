

# Polygon summarization

# this script calculates latitudinal and elevational bounds 
# of individual study polygons


poly_sum <- function(data){
    
    if(nrow(data)>1) stop("this demo runs a single polygon at a time.")

    # fix geometry if needed:
    fc <- sf::st_make_valid(data)
    
    #project Mollweide (less deformation)
    fc <- terra::vect(fc)
    fcMoll <- terra::project(fc, Moll)
    
    
    #summarize attribute information
    area <- sum(terra::expanse(fcMoll, unit = "km"))
    minLat <- terra::ext(fcMoll)[3]
    maxLat <- terra::ext(fcMoll)[4]
    extendLat <- abs(maxLat - minLat)/1000
    centroids <- terra::centroids(fc)
    centroids <- terra::geom(centroids)
    
    

    if( any(data$type == "ELE")){ #if study area linked to elevation shifts
        elevation_tmp <- terra::crop(elevationM, fcMoll)
        Eleminm <- terra::global(elevation_tmp, "min", na.rm = TRUE)[,1]
        Elemaxm <- terra::global(elevation_tmp, "max", na.rm = TRUE)[,1]
        Elemeanm <- terra::global(elevation_tmp, "mean", na.rm = TRUE)[,1]
        extendEle <- Elemaxm - Eleminm
        
    } else{
        Eleminm <- NA 
        Elemaxm <- NA 
        Elemeanm <- NA 
        extendEle <- NA
    }
    
    
    # create df
    met <- tibble( 
        article_id =  data$article_id,
        poly_id = data$poly_id,
        lat_cent_deg = round(centroids[4],5),
        lon_cent_deg = round(centroids[3],5),
        area_km2 = round(as.numeric(area),2),
        lat_extent_km = round(extendLat,2),
        lat_min = round(minLat,2),
        lat_max = round(maxLat,2),
        ele_extent_m = round(extendEle,2),
        ele_min_m = round(Eleminm,2),
        ele_max_m = round(Elemaxm,2),
        ele_mean_m = round(Elemeanm,2),
    )
    
    # add "demo" to colnames
    colnames(met)[3:12] <- paste0(colnames(met)[3:12],"_demo")

    data <- data %>%
        left_join(met)
    
}





