# Interpolate DUALEM measurements

library(terra)
library(magrittr)
library(viridis)
library(tidyr)
library(tidyterra)
library(tools)
library(fields) #tps
library(dplyr)
library(sf)
library(gstat)
library(automap)
library(obliquer)
library(ranger)
library(caret)
library(tibble)

if (!exists("started")) {
  wd <- getwd()
  setwd("..")
  root <- getwd()
  setwd(wd)
  rm(wd)
  started <- TRUE
}

source("f_move_pts_time.R")

dir_dat <- paste0(root, "/PIbalance_gis/")
# sitenames <- c("Field_3", "Field_39", "Field_47", "Trial_070552323", "Vindum")
sitenames <- c("Field_3", "Field_39", "Field_47", "Vindum")

move_values <- c(1, 3, 1, 1)
# When moving points, remove points with large gaps (breaks/outliers).
# Maybe move points in the direction of the next point? Use time diff?

pts_per_gap_width <- c(6, 16, 4, 3)

# The real gap ratio for Field 47 is 3.4

# Field 47: Remove points affected by power lines
# First move points, then remove the outliers
# For GPS time between 132000 and 136000, points should not be moved.

# Do parcels later

sites_df <- data.frame(name = sitenames)

sites_df %<>% mutate(
  filename = "DUALEM1.shp",
  channel = "AUX_X4",
  filename = case_when(
    name == "Vindum" ~ "DUALEM21.shp",
    .default = filename
  ),
  channel = case_when(
    name == "Vindum" ~ "PRP1m",
    .default = channel
  )
)

dir_results <- root %>%
  paste0(., "/results/") %T>%
  dir.create()

# Loop

set.seed(1234)

seeds <- sample(c(10000:99999), 1000)

# for (i in 1:nrow(sites_df)) {
for (i in 2) {
  
  dir_cov <- paste0(dir_dat, sites_df$name[i], "/covariates/")
  
  dem <- dir_cov %>%
    paste0(., "dem_2m.tif") %>%
    rast()
  
  ogcs <- obliquify(dem, 64)
  
  dualem <- paste0(
    dir_dat, sites_df$name[i], "/DUALEM_raw/", sites_df$filename[i]
  ) %>%
    vect()
  
  if (move_values[i] > 0) {
    if (sites_df$name[i] == "Field_47") {
      # Identify points to remove
      
      dualem_crop_out <- paste0(
        dir_dat, sites_df$name[i], "/DUALEM_raw/DUALEN_Crop_out.shp"
      ) %>%
        vect()
      
      dualem$overlap <- is.related(dualem, dualem_crop_out, "within")
      
      # Move points 
      
      dualem_move <- dualem %>%
        filter(GPS_TIME < 130000 | GPS_TIME > 136000)

      dualem_nomove <- dualem %>%
        filter(GPS_TIME > 130000,
               GPS_TIME < 136000)

      dualem_moved1 <- move_pts_time(
        dualem_move,
        time_shift = -0.5*move_values[i]
      )
      
      dualem_moved2 <- move_pts_time(
        dualem_nomove,
        time_shift = 0.5
      )

      dualem_moved <- c(
        dualem_moved1,
        dualem_moved2
        # dualem_nomove
      ) %>%
        vect()
      
      # dualem_moved <- move_pts_time(
      #   dualem,
      #   time_shift = -0.5*move_values[i]
      # )
      
      # Remove points
      
      dualem_moved %<>% filter(overlap == FALSE)
      
    } else {
      if (sites_df$name[i] == "Vindum") {
        dualem_moved <- dualem[-c(1:move_values[i]), ] %>%
          geom() %>%
          vect()
        
        crs(dualem_moved) <- crs(dualem)
        
        values(dualem_moved) <- dualem %>%
          values() %>%
          slice_head(n = nrow(.) - move_values[i])
      } else {
        dualem_moved <- move_pts_time(
          dualem,
          time_shift = -0.5*move_values[i]
        )
      }
    }
  } else {
    dualem_moved <- dualem
  }
  
  # Interpolation using random forest
  
  # pts1 <- terra::extract(ogcs, dualem_moved, bind = TRUE) %>%
  #   values() %>%
  #   select(any_of(c(sites_df$channel[i], names(ogcs)))) %>%
  #   drop_na()
  # 
  # fm1 <- as.formula(
  #   paste0(sites_df$channel[i], "~", paste(names(ogcs), collapse = ' + '))
  # )
  # 
  # bootControl <- trainControl(method = "oob", number = 1, verboseIter = TRUE)
  # 
  # rfGrid <- expand.grid(
  #   mtry = 8,
  #   splitrule = c("extratrees"),
  #   min.node.size = 1
  # )
  # 
  # set.seed(1)
  # 
  # ogcmodel <- train(
  #   fm1, 
  #   pts1,
  #   method = "ranger",
  #   trControl = bootControl
  #   ,
  #   tuneGrid = rfGrid,
  #   sample.fraction = 1/pts_per_gap_width[i],
  #   num.trees = 500*pts_per_gap_width[i]
  # )
  # 
  # ogcmodel
  # 
  # p_ogc <- predict(ogcs, ogcmodel, na.rm = TRUE)
  # 
  # plot(p_ogc)
  # 
  # names(p_ogc) <- "DUALEM_PRP1m"
  # 
  # writeRaster(
  #   p_ogc,
  #   filename = paste0(dir_cov, "DUALEM_PRP1m.tiff"),
  #   overwrite = TRUE
  #   )
  
  
  #################
  
  # General stuff
  
  dem_ext_rast <- dem %>%
    ext() %>%
    rast(
      resolution = res(dem),
      crs = crs(dem)
    )
  
  dualem_moved %<>%
    select(any_of(sites_df$channel[i])) %>% unique()
  
  dualem_df <- dualem_moved %>%
    select(any_of(sites_df$channel[i])) %>%
    as.data.frame(geom = "XY")
  
  # Using kriging
  
  # forced_nugget <- dualem %>%
  #   select(any_of(sites_df$channel[i])) %>%
  #   values() %>%
  #   unlist() %>%
  #   diff() %>%
  #   raise_to_power(2) %>%
  #   mean()
  
  field_area <- expanse(dem, transform = FALSE) %>%
    select(area) %>%
    unlist() %>%
    unname()
  
  n_pts <- terra::extract(dem, dualem_moved) %>%
    na.omit() %>%
    nrow()
  
  dens_expected <- n_pts / field_area
  
  sigma_pts <- sqrt(field_area / (n_pts * pi))
  
  library(spatstat.geom)
  library(spatstat)
  dens_out <- ppp(
    dualem_df[, 2],
    dualem_df[, 3],
    ext(dualem_moved)[1:2],
    ext(dualem_moved)[3:4]
  ) %>%
    density(
      sigma = sigma_pts,
      at = "points",
      leaveoneout = FALSE
    )
  
  attributes(dens_out) <- NULL
  
  weights_dens <- dens_expected / dens_out
  
  weights_dens[weights_dens > 1] <- 1
  
  sum(weights_dens)
  
  plist <- list()
  
  for (j in 1:100) {
    set.seed(seeds[j])
    
    dualem_bootsample <- dualem_moved %>%
      sample(
        round(sum(weights_dens)) / (pts_per_gap_width[i]),
        prob = weights_dens
      )
    
    id_smalldist <- distance(dualem_bootsample) %>%
      apply(1, function(x) {sum(x < (sigma_pts))}) %>%
      unname()
    
    s <- dualem_bootsample %>%
      select(any_of(sites_df$channel[i])) %>%
      sf::st_as_sf(.)
    
    variogram <- autofitVariogram(
      as.formula(paste0(sites_df$channel[i], " ~ 1")),
      s
      # ,
      # fix.values = c(forced_nugget, NA, NA)
    )
    
    xy <- geom(dualem_bootsample) %>%
      as.data.frame() %>%
      select(x, y)
    
    gOK <- gstat(
      NULL,
      sites_df$channel[i],
      as.formula(paste0(sites_df$channel[i], " ~ 1")),
      dualem_bootsample %>%
        select(any_of(sites_df$channel[i])) %>%
        as.data.frame(geom = "XY"),
      locations = ~ x + y,
      model = variogram$var_model,
      # nmin = 6,
      maxdist = sqrt(field_area)/2,
      omax = 8,
      force = TRUE
    )
    
    lower_p <- max(c(0, mean(s[[1]]) - sd(s[[1]])*3))
    upper_p <- mean(s[[1]]) + sd(s[[1]])*3
    
    p <- interpolate(dem_ext_rast, gOK, index = 1)
    p <- mask(p, dem)
    p <- clamp(p, lower = lower_p, upper_p)
    
    plist[[j]] <- p
  }
  pkrig <- plist %>% rast() %>% mean(na.rm = TRUE)
  
  plot(pkrig)
  
  writeRaster(
    pkrig,
    filename = paste0(dir_cov, "DUALEM_PRP1m.tif"),
    names = "DUALEM_PRP1m",
    overwrite = TRUE
    )


  #########
  
  # Fast Tps
  
  # dualem_distinct <- dualem_df %>%
  #   mutate(
  #     x_round = round(x),
  #     y_round = round(y)
  #   ) %>%
  #   distinct(x_round, y_round, .keep_all = TRUE)
  # 
  # set.seed(1)
  # 
  # dualem_bootsample <- dualem_distinct %>%
  #   sample_frac(
  #     size = 1 / (pts_per_gap_width[i] * 2)
  #     )
  # 
  # xy <- dualem_bootsample %>%
  #   select(x, y) %>%
  #   as.matrix()
  # 
  # v <- dualem_bootsample %>%
  #   select(any_of(sites_df$channel[i])) %>%
  #   unlist() %>%
  #   unname()
  # 
  # library(fields)
  # tps <- fastTps(xy, v, aRange = 30)
  # p <- dem_ext_rast
  # 
  # p <- interpolate(p, tps)
  # p <- mask(p, dem)
  # plot(p)
}



# Test moved points and cropping

# plot(dualem, "AUX_X4")
# 
# plot(dualem_moved, "AUX_X4")
# 
# i <- 3
# 
# dualem <- paste0(
#   dir_dat, sites_df$name[i], "/DUALEM_raw/", sites_df$filename[i]
# ) %>%
#   vect()
# 
# dualem_crop_out <- paste0(
#   dir_dat, sites_df$name[i], "/DUALEM_raw/DUALEN_Crop_out.shp"
# ) %>%
#   vect()
# 
# dualem$overlap <- is.related(dualem, dualem_crop_out, "within")
# 
# plot(dualem, "overlap")
# 
# plot(dualem$GPS_TIME)
# 
# plot(dualem$GPS_TIME, dualem$ID2)
# 
# plot(dualem$GPS_TIME)
# 
# xy_distances <- dualem %>%
#   arrange(GPS_TIME) %>%
#   geom() %>%
#   as.data.frame() %>%
#   select(x, y) %>%
#   as.matrix() %>%
#   diff() %>%
#   as.data.frame() %>%
#   mutate(dist = sqrt(x^2 + y^2))
# 
# time_diff <- dualem$GPS_TIME %>% sort() %>% diff()
# 
# mean(time_diff)
# sd(time_diff)
# 
# median(time_diff)
# 
# plot(time_diff, ylim = c(0, 20))
# 
# plot(dualem$GPS_TIME[-1], time_diff, ylim = c(0, 20))
# 
# plot(dualem$GPS_TIME[-1], time_diff, ylim = c(0, 20), xlim = c(129000, 137000))
# 
# data.frame(
#   GPS_TIME = sort(dualem$GPS_TIME[-1]),
#   time_diff = time_diff
#   ) %>%
#   filter(time_diff > 100)
# 
# plot(xy_distances$dist / time_diff)
# 
# speed <- xy_distances$dist / time_diff
# 
# accelleration <- diff(speed)
# 
# plot(accelleration)
# 
# plot(xy_distances$dist)
# 
# plot(xy_distances$dist, ylim = c(0, 20))
# 
# xy_distances$dist %>% mean()
# 
# sd(xy_distances$dist)
# 
# mean(xy_distances$dist) + sd(xy_distances$dist)*3
# 
# xy_distances$dist %>%
#   diff() %>%
#   plot()
# 
# xy_distances$dist %>%
#   diff() %>%
#   mean()
# 
# xy_distances$dist %>%
#   diff() %>%
#   sd()

# Try moving points based on time

pts_try <- dualem[1:90] %>%
  arrange(GPS_TIME)

plot(pts_try)

source("f_move_pts_time.R")

i <- 1

dualem <- paste0(
  dir_dat, sites_df$name[i], "/DUALEM_raw/", sites_df$filename[i]
) %>%
  vect()

xy_distances <- dualem %>%
  arrange(GPS_TIME) %>%
  geom() %>%
  as.data.frame() %>%
  select(x, y) %>%
  as.matrix() %>%
  diff() %>%
  as.data.frame() %>%
  mutate(dist = sqrt(x^2 + y^2))

mean(xy_distances$dist)
sd(xy_distances$dist)

mean(xy_distances$dist) + sd(xy_distances$dist)*3

time_diff <- dualem$GPS_TIME %>% sort() %>% diff()

mean(time_diff)
sd(time_diff)

median(time_diff)

pts_try_moved <- move_pts_time(
  dualem,
  time_shift = -0.5
)

pts_try_moved

values(pts_try_moved)

plot(pts_try_moved$xy_gap)

plot(dualem, "AUX_X4")
plot(pts_try_moved, "AUX_X4")

pts_try_moved$dist_moved %>% plot()

# For small plot
# plot(pts_try)
# plot(pts_moved_vect, add = TRUE, col = "red", pch = 1) 

## Stuff for kriging


# dualem_df <- dualem %>%
#   select(any_of(sites_df$channel[i])) %>%
#   as.data.frame(geom = "XY")
# 
# dem_ext_rast <- dem %>%
#   ext() %>%
#   rast(
#     resolution = res(dem),
#     crs = crs(dem)
#   )

# forced_nugget <- dualem %>%
#   select(any_of(sites_df$channel[i])) %>%
#   values() %>%
#   unlist() %>%
#   diff() %>%
#   raise_to_power(2) %>%
#   mean()
# 
# s <- dualem %>%
#   select(any_of(sites_df$channel[i])) %>%
#   sf::st_as_sf(.)
# 
# variogram <- autofitVariogram(
#   as.formula(paste0(sites_df$channel[i], " ~ 1")),
#   s,
#   fix.values = c(forced_nugget, NA, NA)
#   )
# 
# xy <- geom(dualem) %>%
#   as.data.frame() %>%
#   select(x, y)
# 
# gOK <- gstat(
#   NULL,
#   sites_df$channel[i],
#   as.formula(paste0(sites_df$channel[i], " ~ 1")),
#   dualem %>%
#     select(any_of(sites_df$channel[i])) %>%
#     as.data.frame(geom = "XY"),
#   locations = ~ x + y,
#   model = variogram$var_model,
#   nmin = 8,
#   maxdist = 30,
#   omax = 4,
#   force = TRUE
#   )
# 
# p <- interpolate(dem_ext_rast, gOK, index = 1)
# p <- mask(p, dem)
# plot(p)

## Fast Tps

# dualem_distinct <- dualem_df %>%
#   mutate(
#     x_round = round(x),
#     y_round = round(y)
#   ) %>%
#   distinct(x_round, y_round, .keep_all = TRUE)
# 
# xy <- dualem_distinct %>%
#   select(x, y) %>%
#   as.matrix()
# 
# v <- dualem_distinct %>%
#   select(any_of(sites_df$channel[i])) %>%
#   unlist() %>%
#   unname()
# 
# library(fields) 
# tps <- fastTps(xy, v, aRange = 30)
# p <- dem_ext_rast
# 
# p <- interpolate(p, tps)
# p <- mask(p, dem)
# plot(p)

# END