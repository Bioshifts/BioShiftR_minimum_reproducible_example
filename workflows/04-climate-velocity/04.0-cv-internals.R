
# climate velocity -- necessary functions

# this script calculates latitudinal or elevational climate velocity from
# climate raster data. 

# Note that these functions were based on the functions from:
# 1. the climetrics package (https://github.com/shirintaheri/climetrics)
# 2. gVelocity workflow (functions `tempTrend`, `spatGrad`, `gVoCC`) from 
#    the (depreciated) VoCC package (https://github.com/JorGarMol/VoCC)
# which we adapted to achieve poleward standardization of gradient-based cvs.


# temp trend internal helpers ---------------------------------------------
# extract slope coefficient of linear model
# inherited from climetrics:
temp_gradFun <- function(x, th) {
    x <- x[!is.na(x)]
    if (length(x) > th) {
        s <- lm(x~c(1:length(x)))
        s$coefficients[2] 
    } else NA
}

# temp trend function -----------------------------------------------------
# function to apply temp_gradFun to every raster cell
# inherited from climetrics:
temp_trend <- function(x, th, ncores=NULL, file_name="",overwrite=TRUE) {
    if(is.null(ncores)){
        tmp <- terra::app(x,temp_gradFun,th=th,filename=file_name,overwrite=overwrite)
    } else {
        tmp <- terra::app(x,temp_gradFun,th=th,cores=ncores,filename=file_name,overwrite=overwrite)
    }
    names(tmp) <- "Trend" # name the layer
    return(tmp)
}



# spatial gradient internal helpers ---------------------------------------

#----
# Utils
# inherited from climetrics:
.is_package_installed <- function(n) {
    names(n) <- n
    sapply(n, function(x) length(unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))) > 0)
}
#----
# inherited from climetrics:
# get projection type
.getProj <- function(x) {
    if (inherits(x,'Raster')) {
        if (!is.na(projection(x))) strsplit(strsplit(projection(x),'\\+proj=')[[1]][2],' ')[[1]][1]
        else {
            if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
            else 'projected'
        }
    } else {
        if (!is.na(crs(x))) strsplit(strsplit(crs(x,proj=TRUE),'\\+proj=')[[1]][2],' ')[[1]][1]
        else {
            if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
            else 'projected'
        }
    }
}
#-----
# inherited from climetrics:
# weighted mean of 6 pairwise comparisons in a 9-cell neighborhood
.mnwm <- function(d1, d2, d3, d4, d5, d6){
    X <- sum(c(d1, d2*2, d3, d4, d5*2, d6), na.rm = T)
    w <- sum(c(1,2,1,1,2,1) * is.finite(c(d1, d2, d3, d4, d5, d6)))
    return(X/w)
}
#-----
# inherited from climetrics:
# find angle of warmer temperature
.ang <- function(dx, dy){
    ifelse(dy < 0, 180 + rad_to_deg(atan(dx/dy)),
           ifelse(dx < 0, 360 + rad_to_deg(atan(dx /dy )), rad_to_deg(atan(dx/dy))))
}


#---
# change from degrees to radians
deg_to_rad <- function (degree) {
    (degree * pi) / 180
}
#---
# change from radians to degrees
rad_to_deg <-  function (radian) {
    (radian * 180) / pi
}
#---
d2km <- function (d, base.latitude = 1) 
{
    if (!requireNamespace("fields")) 
        stop("Required fields package is missing.")
    onerad_to_degree.dist <- fields::rdist.earth(matrix(c(0, base.latitude), ncol = 2), 
                                                 matrix(c(1, base.latitude), ncol = 2), 
                                                 miles = FALSE)[,1]
    out <- d * onerad_to_degree.dist
    return(out)
}



#--------------
# spatial gradient
# function to find the spatial gradient of mean temperature across space.
# inherited from climetrics with small adaptation to include neighboring cells
# in the polygon buffer when calculating spatial gradient for focal cells, but
# then remove them from the final raster
spatial_grad <- function(rx, 
                         y_diff = 1, 
                         unit_out = "km") {
    
    if(!(unit_out == "km" | unit_out == "m")){
        stop("unit_out should be 'km' or 'm'")
    }
    
    # if there's more than one layer (i.e., input data weren't averaged), average them.
    if(nlyr(rx) > 1){ rx <- mean(rx,na.rm = TRUE) }
    
    # find y-distance between cells
    if (.getProj(rx) == 'longlat') {
        if(unit_out=="km"){
            y_dist <- d2km(res(rx)) # from degrees to km
        }
        if(unit_out=="m"){
            y_dist <- d2km(res(rx))*1000 # from degrees to m
        }
        
    } else {
        if(unit_out=="km"){
            y_dist <- res(rx) / 1000 # from meters to km
        }
        if(unit_out=="m"){
            y_dist <- res(rx)
        }
        y_diff <- NA
    }
    
    if (!.is_package_installed("dplyr") || !.is_package_installed('tidyr')) stop('The packages dplyr and tidyr are needed for this metric; Please make sure they are installed!')
    
    # find cells that touch each other
    y <- data.frame(adjacent(rx, cells=1:ncell(rx), directions=8,pairs=TRUE))
    y <- y[order(y$from, y$to),] # sort 
    y <- na.omit(y) # remove NAs (but there shouldn't be any)
    
    
    # add mean temp from each "to" cell
    y$temp <- rx[y$to][,1]
    # find the difference in rows between adjacent cells
    # -1 = focal cell is above; 0 = same row; 1 = focal cell is below
    y$sy <- rowFromCell(rx, y$from)-rowFromCell(rx, y$to)
    y$sx <- colFromCell(rx, y$to)-colFromCell(rx, y$from)
    y$sx[y$sx > 1] <- -1
    y$sx[y$sx < -1] <- 1
    
    # paste the coordinates of the col and row difference
    y$code <- paste(y$sx, y$sy)
    
    # translate coordinates to plain English
    y$code1 <- eval(parse(text='dplyr::recode(y$code,
                           `1 0` = "tempE",
                           `-1 0` = "tempW",
                           `-1 1` = "tempNW",
                           `-1 -1` = "tempSW",
                           `1 1` = "tempNE",
                           `1 -1` = "tempSE",
                           `0 1` = "tempN",
                           `0 -1` = "tempS")'),envir =environment())
    
    # remove +1, 0, -1 directional columns, keep only plain english. 
    y3b <- eval(parse(text="dplyr::select(y,from, code1, temp)"),envir =environment())
    # spread so each focal cell is one row
    y3b <- eval(parse(text="tidyr::spread(y3b,code1, temp)"),envir =environment())
    # add temperature of the "from" cells. 
    y3b$tempFocal <- rx[y3b$from][,1]
    # add latitude (used to correct for differences in cell area across latitudes.)
    y3b$LAT <- yFromCell(rx, y3b$from)
    
    if(!is.na(y_diff)) {
        y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = cos(deg_to_rad(LAT + y_diff)),
                         latneg = cos(deg_to_rad(LAT - y_diff)),
                         latfocal = cos(deg_to_rad(LAT)))"),envir =environment())
    } else {
        
        y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = 1,
                         latneg = 1,
                         latfocal = 1)"),envir =environment())
    }
    
    # calculate the difference in temperature between 9 neighboring cells.
    y3c <- "dplyr::mutate(y3b,
                       gradWE1 = (tempN-tempNW)/ (latpos *  y_dist[1]),
                       gradWE2 = (tempFocal - tempW)/(latfocal * y_dist[1]),
                       gradWE3 = (tempS-tempSW)/(latneg * y_dist[1]),
                       gradWE4 = (tempNE-tempN)/(latpos * y_dist[1]),
                       gradWE5 = (tempE-tempFocal)/(latfocal * y_dist[1]),
                       gradWE6 = (tempSE-tempS)/(latneg * y_dist[1]),
                       gradNS1 = (tempNW-tempW)/y_dist[2],
                       gradNS2 = (tempN-tempFocal)/y_dist[2],
                       gradNS3 = (tempNE-tempE)/y_dist[2],
                       gradNS4 = (tempW-tempSW)/y_dist[2],
                       gradNS5 = (tempFocal-tempS)/y_dist[2],
                       gradNS6 = (tempE-tempSE)/y_dist[2])" 
    
    y3c <- eval(parse(text=y3c),envir=environment())
    
    y3c <- eval(parse(text="dplyr::rowwise(y3c)"),envir=environment())
    
    y3c <- eval(parse(text="dplyr::mutate(y3c,
      WEgrad = .mnwm(gradWE1, gradWE2, gradWE3, gradWE4, gradWE5, gradWE6),
      NSgrad = .mnwm(gradNS1, gradNS2, gradNS3, gradNS4, gradNS5, gradNS6),
      angle = .ang(WEgrad, NSgrad))"),envir=environment())
    
    y3c <- eval(parse(text="dplyr::select(y3c,icell = from, WE = WEgrad, NS = NSgrad, angle = angle)"),envir=environment())
    
    # edit from climetrics: replace grad and angle with NA if input raster had no value 
    y3c$angle[is.na(values(rx))] <- NA
    y3c$NS[is.na(values(rx))] <- NA
    y3c$WE[is.na(values(rx))]  <- NA
    
    NS <- y3c$NS # isolate NS gradient 
    WE <- y3c$WE # isolate EW gradient
    NS[is.na(NS)] <- 0 # change NAs to 0
    WE[is.na(WE)] <- 0 # change NAs to 0
    NAsort <- ifelse((abs(NS)+abs(WE)) == 0, NA, 1)
    # calculate spatial gradient (hypotneuse of NS and EW gradient)
    y3c$Grad <- NAsort * sqrt((WE^2) + (NS^2))
    
    
    # make output raster
    rx <- c(rx, rx, rx, rx)
    names(rx) <- names(y3c)[-1]
    rx[[1]][y3c$icell] <- y3c$WE
    rx[[2]][y3c$icell] <- y3c$NS
    rx[[3]][y3c$icell] <- y3c$angle
    rx[[4]][y3c$icell] <- y3c$Grad
    
    return(rx)
}


# calculate gradient-based climate velocity.
# adapted from climetrics to control for cooling cells. Because gradient is 
# always positive (calculated as a hypotneuse), 
gVelocity <- function(grad, # spatial gradient raster
                         slope, # temp trend raster
                         grad_col = "Grad", 
                         truncate=FALSE) {
    
    v <- slope # this is temp trend 
    g <- grad[grad_col] # this is gradient == always positive
    v_ang <- grad["angle"] # angle points to direction of WARMING
    
    # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
    v_ang <- terra::ifel(slope > 0, v_ang + 180, v_ang) %% 360
    # NOTE: this line is inherited from VoCC/climetrics, but we make some changes below. 
    # this above comment is true: gVel angles should be opposite to the spatial gradient 
    # if the trend is warming and same angle as spatial gradient if cooling, 
    # but this strategy omits that the numerator of CV (ttrend) can be negative, 
    # and thus influences the interpretation of the ending vector. 
    
    # For example: --- --- --- ---
    # 1. a cell with a southward gradient (angle = 180; pointing towards warmer south), 
    #    and a temperature trend of 1 (e.g., warming) will result in a vector with 
    #    magnitude 1 and angle 0 (towards North). That is correct. 
    # 2. However, a cell with a southward gradient (angle = 180; pointing towards warmer south) 
    #    and a temperature trend of -1 (e.g., cooling) *should* result in a vector of 
    #    magnitude 1 and angle 180 (south) (to represent southward movement of isotherms),
    #    but the current code results in a magnitude -1 and angle 180, which is 
    #    equivalent to a northward velocity (mag 1 at angle 0), which does not 
    #    accurately represent the movement of isotherms, in this case. 
    
    # Potential fixes: --- --- --- ---
    # 1. keep current code (flipping angle when trend is positive, as done above), 
    #    but then use abs(v)/g to calculate the vector magnitude
    #    this way, velocities are always positive, and angles always point in the 
    #    direction of isotherm movement. 
    # 2. standardize all spgrad values to direction of cooling by flipping ALL 
    #    angles, not just when ttrend is positive.
    
    # Becasue since we already flipped half the angles in the line above (inherited),
    # we'll implement Fix 1 here and diverge from VoCC/climetrics.

    # calculate velocity
    # note that this differs from VoCC/climetrics
    v <- abs(v) / g
    
    # truncate extreme values if truncate=T 
    # (inherited from climetrics/VoCC)
    if(truncate){
        .o <- as.matrix(global(v,fun=quantile,probs=c(0.05,0.95),na.rm=TRUE))[1,]
        v[v < .o[1]] <- .o[1]
        v[v > .o[2]] <- .o[2] 
    }
    # note that anomalously high climate velocity values often occur when cells in a 
    # 9-cell neighborhood used to calculate spatial gradient are NA. Thus, our
    # addition of a buffer around the study polygon helps to reduce anomalous values
    # within the focal polygon by keeping buffer cells, and making them "sacrificial"
    # (e.g., they might have anomalously high values due to their surrounding NAs,
    # but won't - or will barely - be averaged in to the study area mean value; and
    # their inclusion in the spatial gradient calculation will mean that focal cells
    # within the study area less often have surrounding NAs in their neighborhoods).

    # export velocity magnitude and velocity angle layers
    output <- c(v,v_ang)
    names(output) <- c("Vel", "Ang")
    return(output)
}




#----------
# version for FLIPPED spatial gradients, so we can omit the absolute value and angle flipping here.
# truncate is for bounding max and min values to upper (95%) and lower (5%) quantiles, respectively
gVelocity_JL_flipped <- function(grad, slope, grad_col = "Grad", truncate=FALSE) {
    
    v <- slope # this is temp trend 
    g <- grad[grad_col] # this is gradient == always positive
    v_ang <- grad["angle"] # angle points to direction of COOLING
    
    # calculate velocity
    v <- (v / g)
    
    # truncate extreme values if truncate=T (if there are lots of cells in the raster)
    if(truncate){
        .o <- as.matrix(global(v,fun=quantile,probs=c(0.05,0.95),na.rm=TRUE))[1,]
        v[v < .o[1]] <- .o[1]
        v[v > .o[2]] <- .o[2] 
    }

    # # this produces some negative velocities, which are the same as positive velocities
    # # in the opposite direction, so here, flip the angles of the negative velocities 
    # # and take absolute value of the velocity
    v_ang[v < 0] <- (v_ang[v < 0] + 180) %% 360
    v <- abs(v)
    
    # export velocity magnitude and velocity angle layers
    output <- c(v,v_ang)
    names(output) <- c("Vel", "Ang")
    return(output)
}




