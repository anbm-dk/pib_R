# Resample 10 texture maps

library(terra)
library(magrittr)
library(viridis)

if (!exists("started")) {
  wd <- getwd()
  setwd("..")
  root <- getwd()
  setwd(wd)
  rm(wd)
  started <- TRUE
}

dir_dat <- paste0(root, "/PIbalance_gis/")
sitenames <- c("Field_3", "Field_39", "Field_47", "Trial_070552323", "Vindum")

bigmaps <- c(
  paste0(
    "O:/Tech_AGRO/Jord/anbm/GIS/Texture_maps_10m/depth_000_030_cm/",
    "Tekstur2024_000_030_cm/",
    "Clay_mean.tif"
  ),
  paste0(
    "O:/Tech_AGRO/Jord/anbm/GIS/Texture_maps_10m/depth_000_030_cm/",
    "Tekstur2024_000_030_cm/",
    "Silt_mean.tif"
  ),
  paste0(
    "O:/Tech_AGRO/Jord/anbm/GIS/Texture_maps_10m/depth_000_030_cm/",
    "Tekstur2024_000_030_cm/",
    "Fine_sand_mean.tif"
  ),
  paste0(
    "O:/Tech_AGRO/Jord/anbm/GIS/Texture_maps_10m/depth_000_030_cm/",
    "Tekstur2024_000_030_cm/",
    "Coarse_sand_mean.tif"
  ),
  paste0(
    "O:/Tech_AGRO/Jord/anbm/GIS/Texture_maps_10m/depth_000_030_cm/",
    "Kulstof2022_000_030_cm/Kulstof_kombineret.tif"
  )
)

outnames <- c("clay", "silt", "fine_sand", "coarse_sand", "SOC") %>%
  paste0(., "_10m")


for (i in 1:length(sitenames)) {
  
  outmaps_i <- list()
  
  for (j in 1:length(bigmaps)) {
    
    dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")
    
    dem_square <- paste0(dir_dat, sitenames[i], "/dem_2m/dem_2m.tif") %>%
      rast()
    
    bigmap <- bigmaps[j] %>%
      rast()
    bigmap_crop <- crop(bigmap, dem_square)
    
    my_focal_weights <- focalMat(
      bigmap_crop,
      c(10, 20),
      type = c('Gauss')
    )
    
    r1 <- focal(bigmap_crop, w = my_focal_weights, na.rm = TRUE)
    r2 <- focal(is.finite(bigmap_crop), w = my_focal_weights, na.rm = TRUE)
    r3 <- r1/r2
    
    dem_r <- dir_cov %>% 
      paste0(., "/dem_2m.tif") %>%
      rast()
    
    r_res <- resample(r3, dem_r, method = "cubicspline") %>%
      mask(dem_r)
    
    outfile <- dir_cov %>%
      paste0(., outnames[j], ".tif")
    
    names(r_res) <- outnames[j]
     
    writeRaster(
      r_res, filename = outfile, overwrite = TRUE
    )
    
    outmaps_i[[j]] <- r_res
  }
  
  outmaps_i %<>% rast()
  
  names(outmaps_i) <- outnames
  
  plot(outmaps_i, col = cividis(100))
  
  plot(dem_r, col = cividis(100))
}

# Red band for Vindum

dir_cov <- paste0(dir_dat, "Vindum", "/covariates/")

dem_vindum <- dir_cov %>% 
  paste0(., "/dem_2m.tif") %>%
  rast()

dem_vindum_fine <- disagg(dem_vindum, 10)

red_files <- dir_dat %>%
  paste0(., "/Vindum/ortho_red_raw/") %>%
  paste0(., c("Vindum_red_raw1.tif", "Vindum_red_raw2.tif"))
  
red_outmaps <- list()

for (i in 1:length(red_files)) {
  r <- red_files[i] %>% rast()
  
  r_resample <- r %>% resample(dem_vindum_fine)
  
  crs(r_resample) <- crs(dem_vindum)
  
  r_2m <- r_resample %>% aggregate(10) %>% crop(dem_vindum, mask = TRUE)
  
  red_outmaps[[i]] <- r_2m
}

outname <- "ortho_red"
outfile <- dir_cov %>% paste0(., outname, ".tif")

red_out <- red_outmaps %>% rast() %>% mean() %>% round(digits = 0)

names(red_out) <- outname

writeRaster(
  red_out, filename = outfile, overwrite = TRUE
)

# xy rasters

library(obliquer)

outnames <- c("UTMX", "UTMY")

for (i in 1:length(sitenames)) {
  
  dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")
  
  dem_r <- dir_cov %>% 
    paste0(., "/dem_2m.tif") %>%
    rast()
  
  xy_r <- obliquify(dem_r, n_angles = 2)
  
  
  names(xy_r) <- outnames
  
  for (j in 1:length(outnames)) {
    outfile <- dir_cov %>%
      paste0(., outnames[j], ".tif")
    
    writeRaster(xy_r[[j]], filename = outfile, overwrite = TRUE)
  }
  
  plot(xy_r)
}


# END