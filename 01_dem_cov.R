# Generating DEM derivatives for covariates

library(terra)
library(magrittr)
library(viridis)
library(RSAGA)
library(dplyr)

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

for (i in 1:length(sitenames)) {
  dem_r <- paste0(dir_dat, sitenames[i], "/dem_2m/dem_2m.tif") %>%
    rast()
  
  field_polygon <- paste0(dir_dat, sitenames[i], "/field_polygon.shp") %>%
    vect()
  
  obs <- paste0(dir_dat, sitenames[i], "/observations.gpkg") %>%
    vect()
  
  plot(dem_r, col = cividis(100))
  polys(field_polygon)
  plot(obs, "P", add = TRUE, plg = list(x = "topright", bty = "o"))
  
  # Try RSAGA modules
  
  dir <- paste0(dir_dat, sitenames[i], "/dem_2m/")
  DEM <- "dem_2m"
  
  dem_t <- paste0(dir, "/dem_2m.tif") %>%
    rast()
  
  rsaga.env(workspace = dir)
  
  work_env <- rsaga.env(workspace = dir)
  
  rsaga.import.gdal(paste0(DEM, ".tif"), env = work_env)
  
  # Compund analysis - module 0
  
  rsaga.geoprocessor(
    "ta_compound",
    module = 0,
    env = work_env,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      SINKS = "sinks_raw.sgrd", # Depth of sinks
      VALL_DEPTH = "depth_valley.sgrd", # Valley Depth Grid
      RSP = "position_slope_relative.sgrd" # Relative slope position
    )
  )
  
  # Calculate channel depth with dedicated module to keep it below the surface
  
  # Fill no data cells in sinks raster
  
  rsaga.geoprocessor(
    "grid_tools",
    module = 15,
    env = work_env,
    param = list(
      INPUT = "sinks_raw.sgrd",
      RESULT = "sinks.sgrd",
      NODATAOPT = 1,
      NODATA = 0,
      RESULT_NODATA_CHOICE = 1,
      RESULT_NODATA_VALUE = -9
    )
  )
  
  convertibles <- c("sinks", "depth_valley", "position_slope_relative")
  
  for (j in 1:length(convertibles)) {
    rsaga.geoprocessor(
      "io_gdal",
      module = 2,
      env = work_env,
      param = list(
        GRIDS = list(paste0(convertibles[j], ".sgrd")),
        FILE = paste0(convertibles[j], ".tif")
      )
    )
  }
  
  # A: MORPHOMETRY
  
  # Module 0: Slope, Aspect, Curvature
  
  outnames_0 <- c(
    "slope",
    "curvature_general",
    "curvature_profile",
    "curvature_plan",
    "curvature_tangential",
    "curvature_longitudinal",
    "curvature_cross_sectional",
    "curvature_minimal",
    "curvature_maximal",
    "curvature_total",
    "curvature_flow_line"
  )
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 0,
    env = work_env,
    cores = 11,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      SLOPE = paste0(outnames_0[[1]], ".sgrd"),
      C_GENE = paste0(outnames_0[[2]], ".sgrd"),
      C_PROF = paste0(outnames_0[[3]], ".sgrd"),
      C_PLAN = paste0(outnames_0[[4]], ".sgrd"),
      C_TANG = paste0(outnames_0[[5]], ".sgrd"),
      C_LONG = paste0(outnames_0[[6]], ".sgrd"),
      C_CROS = paste0(outnames_0[[7]], ".sgrd"),
      C_MINI = paste0(outnames_0[[8]], ".sgrd"),
      C_MAXI = paste0(outnames_0[[9]], ".sgrd"),
      C_TOTA = paste0(outnames_0[[10]], ".sgrd"), # This one is just 0
      C_ROTO = paste0(outnames_0[[11]], ".sgrd"), # This one is also just 0
      ASPECT = "aspect.sdat"
    )
  )
  
  outr_0 <- outnames_0 %>%
    sapply(function(x) {
      x %>%
        paste0(dir, ., ".sdat") %>%
        rast()
    }) %>%
    rast()
  
  crs(outr_0)
  crs(outr_0) <- crs(dem_t)
  setMinMax(outr_0)
  minmax(outr_0)
  outr_0
  
  outr_0 %>%
    math(.,
         "round",
         digits = 3
    ) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = paste0(dir, outnames_0, ".tif"),
      overwrite = TRUE
    )
  
  # Sine and cosine of aspect
  
  asp_t <- paste0(dir, "aspect.sdat") %>% rast()
  
  crs(asp_t)
  crs(asp_t) <- crs(dem_t)
  setMinMax(asp_t)
  minmax(asp_t)
  asp_t
  
  asp_t %>%
    sin() %>%
    math(.,
         "round",
         digits = 3
    ) %>%
    ifel(is.na(.), 0, .) %>%
    mask(mask = dem_t) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = paste0(dir, "aspect_sine.tif"),
      overwrite = TRUE
    )
  
  asp_t %>%
    cos() %>%
    math(.,
         "round",
         digits = 3
    ) %>%
    ifel(is.na(.), 0, .) %>%
    mask(mask = dem_t) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = paste0(dir, "aspect_cosine.tif"),
      overwrite = TRUE
    )
  
  # OK
  
  # Module 1: Convergence index
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 1,
    env = work_env,
    cores = 11,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      RESULT = "convergence.sgrd"
    )
  )
  
  outr_1 <- paste0(dir, "convergence.sdat") %>% rast()
  
  crs(outr_1)
  crs(outr_1) <- crs(dem_t)
  setMinMax(outr_1)
  minmax(outr_1)
  outr_1
  
  outr_1 %>%
    math(.,
         "round",
         digits = 1,
         filename = paste0(dir, "convergence.tif"),
         overwrite = TRUE,
         datatype = "FLT4S"
    )
  
  # OK
  
  # Modudule 2: Convergence index with search radius
  
  outname_2 <- "convergence_radius"
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 2,
    env = work_env,
    cores = 11,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      CONVERGENCE = "convergence_radius.sgrd"
    )
  )
  
  outr_2 <- paste0(dir, "convergence_radius.sdat") %>% rast()
  
  crs(outr_2)
  crs(outr_2) <- crs(dem_t)
  setMinMax(outr_2)
  minmax(outr_2)
  outr_2
  
  outr_2 %>%
    math(.,
         "round",
         digits = 1,
         filename = paste0(dir, "convergence_radius.tif"),
         overwrite = TRUE,
         datatype = "FLT4S"
    )
  
  # OK
  
  # Module 3: Surface Specific Points (integer format - not ideal)
  # Module 4: Curvature Classification (classes - not ideal)
  # Module 5: Hypsometry (probably not appropriate)
  # Moduel 6: Real Surface Area (probably not useful)
  # Module 7: Morphometric Protection Index (all NA)
  
  # Module 8: Multiresolution Index of Valley Bottom Flatness (MRVBF)
  
  outname_8a <- "MRVBF"
  outname_8b <- "MRRTF"
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 8,
    env = work_env,
    cores = 11,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      MRVBF = paste0(outname_8a, ".sgrd"),
      MRRTF = paste0(outname_8b, ".sgrd")
    )
  )
  
  outr_8a <- outname_8a %>%
    paste0(dir, ., ".sdat") %>%
    rast()
  outr_8b <- outname_8b %>%
    paste0(dir, ., ".sdat") %>%
    rast()
  
  outr_8 <- c(outr_8a, outr_8b)
  
  crs(outr_8)
  crs(outr_8) <- crs(dem_t)
  setMinMax(outr_8)
  minmax(outr_8)
  outr_8
  
  outnames_8 <- outr_8 %>%
    names() %>%
    paste0(dir, ., ".tif")
  
  outr_8 %>%
    math(.,
         "round",
         digits = 2
    ) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = outnames_8,
      overwrite = TRUE
    )
  
  # OK
  
  # Module 9: Downslope Distance Gradient (buggy/unrealistic)
  # Module 10: Mass Balance Index (probably not appropriate)
  
  # Module 11: Effective Air Flow Heights
  # Boehner, J., Antonic, O. (2009): Land-surface parameters specific to
  # topo-climatology. In: Hengl, T., Reuter, H. [Eds.]: Geomorphometry - Concepts,
  # Software, Applications. Developments in Soil Science, Volume 33, p.195-226,
  # Elsevier.
  # https://doi.org/10.1016/S0166-2481(08)00008-1
  
  # This parameter is supposed to act as a proxy for precipitation due to
  # topographic uplift. Slopes facing the wind will receive more precipitation,
  # whereas leewards slopes receive less.
  
  outname_11 <- "height_air_flow"
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 11,
    env = work_env,
    cores = 11,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      AFH = "height_air_flow.sgrd",
      DIR_CONST = 90,
      OLDVER = 1
    )
  )
  
  outr_11 <- paste0(dir, "height_air_flow.sdat") %>% rast()
  
  crs(outr_11)
  crs(outr_11) <- crs(dem_t)
  setMinMax(outr_11)
  minmax(outr_11)
  outr_11
  
  # Apply min filter to remove wierd high outliers (omitted)
  #
  # fmat_min <- matrix(c(1, 1, 1, 1, NA, 1, 1, 1, 1), nrow = 3)
  #
  # r <- outr_11
  #
  # r_mean <- r %>%
  #   focal(
  #     w = fmat_min,
  #     fun = "mean",
  #     na.rm = TRUE,
  #     na.policy = "omit"
  #   )
  #
  # outliers <- numeric()
  #
  # i <- 1
  #
  # outliers[i] <- ifel((r / r_mean) > 3, 1, 0) %>%
  #   values() %>%
  #   sum(na.rm = TRUE) %T>%
  #   print()
  #
  # while (outliers[i] > 0 & i < 1000) {
  #   i %<>% +(1)
  #
  #   r <- ifel(
  #     r_mean != 0,
  #     ifel(
  #       (r / r_mean) > 3,
  #       ifel(
  #         (r - r_mean) > 1,
  #         r_mean,
  #         r
  #       ),
  #       r
  #     ),
  #     0
  #   )
  #
  #   r_mean <- r %>%
  #     focal(
  #       w = fmat_min,
  #       fun = "mean",
  #       na.rm = TRUE,
  #       na.policy = "omit"
  #     )
  #
  #   outliers[i] <- ifel(
  #     r_mean != 0,
  #     ifel(
  #       (r / r_mean) > 3,
  #       ifel(
  #         (r - r_mean) > 1,
  #         1,
  #         0
  #       ),
  #       0
  #     ),
  #     0
  #   ) %>%
  #     values() %>%
  #     sum(na.rm = TRUE) %T>%
  #     print()
  # }
  #
  # outr_11 <- r
  #
  # setMinMax(outr_11, force = TRUE)
  # minmax(outr_11)
  # names(outr_11) <- outname_11
  # outr_11
  #
  # outr_11 <- ifel(
  #   is.na(outr_11),
  #   ifel(
  #     !is.na(dem_t),
  #     0,
  #     NA
  #   ),
  #   outr_11
  # )
  
  outr_11 %>%
    math(.,
         "round",
         digits = 2,
         filename = paste0(dir, "height_air_flow.tif"),
         overwrite = TRUE,
         datatype = "FLT4S"
    )
  
  # OK
  
  # Module 12: Diurnal Anisotropic Heat (probably not useful/appropriate)
  # Module 13: Land Surface Temperature (probably not useful/appropriate)
  
  # Module 14: Relative Heights and Slope Positions
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 14,
    env = work_env
    # , invisible = FALSE
    # , intern = FALSE
    , cores = 11,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      HO = "height_slope.sgrd", # Slope Height
      NH = "height_normalized.sgrd", # Normalized height
      SH = "height_standardized.sgrd", # Standardized Height
      MS = "position_slope_mid.sgrd" # Mid slope position
    )
  )
  
  convertibles <- c(
    "height_slope", "height_normalized", "height_standardized",
    "position_slope_mid"
  )
  
  for (j in 1:length(convertibles)) {
    rsaga.geoprocessor(
      "io_gdal",
      module = 2,
      env = work_env,
      param = list(
        GRIDS = list(paste0(convertibles[j], ".sgrd")),
        FILE = paste0(convertibles[j], ".tif")
      )
    )
  }
  
  # OK (round off results?)
  
  # Module 15: Wind Effect (Windward / Leeward Index)
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 15,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      EFFECT = "wind_effect.sgrd",
      DIR_CONST = 90,
      OLDVER = 1
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "wind_effect.sgrd",
      FILE = "wind_effect.tif"
    )
  )
  
  # OK (round off results´?)
  
  # Module 16: Terrain ruggedness index
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 16,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      TRI = "ruggedness.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "ruggedness.sgrd",
      FILE = "ruggedness.tif"
    )
  )
  
  # Ok (round off results?)
  
  # Module 17: Vector Ruggedness Measure (VRM) [probably not useful]
  
  # Module 18: Topographic Position Index (TPI)
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 18,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      TPI = "position_topographic.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "position_topographic.sgrd",
      FILE = "position_topographic.tif"
    )
  )
  
  # OK (round off results?)
  
  # Module 19: TPI Based Landform Classification (classes not useful)
  # Module 20: Terrain Surface Texture (probably not useful)
  
  # Module 21: Terrain Surface Convexity
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 21,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      CONVEXITY = "convexity.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "convexity.sgrd",
      FILE = "convexity.tif"
    )
  )
  
  # OK (round off decimals?)
  
  # Module 22: Terrain Surface Classification (classes not useful)
  # Module 23: Morphometric Features (redundant)
  # Module 24: Valley and Ridge Detection (Top Hat Approach) (very slow)
  # Module 25: Fuzzy Landform Element Classification (too noisy)
  # Module 26: Upslope and Downslope Curvature (redundant)
  
  # Module 27; Wind Exposition Index
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 27,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      EXPOSITION = "wind_exposition.sgrd",
      STEP = 90,
      OLDVER = 1
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "wind_exposition.sgrd",
      FILE = "wind_exposition.tif"
    )
  )
  
  # OK (round off results?)
  
  # Module 28: Multi-Scale Topographic Position Index (TPI)
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 28,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      TPI = "position_topographic_multi.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "position_topographic_multi.sgrd",
      FILE = "position_topographic_multi.tif"
    )
  )
  
  # OK (round off results?)
  
  # Module 29: Wind Shelter Index
  
  rsaga.geoprocessor(
    "ta_morphometry",
    module = 29,
    env = work_env,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      SHELTER = "wind_shelter.sgrd",
      DIRECTION = 90
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "wind_shelter.sgrd",
      FILE = "wind_shelter.tif"
    )
  )
  
  # OK (round off results?)
  
  # HYDROLOGY
  
  # Module 0: Flow accumulation
  
  rsaga.geoprocessor(
    "ta_hydrology",
    module = 0,
    env = work_env,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      FLOW = "flow_accumulation.sgrd",
      METHOD = 5,
      FLOW_UNIT = 1,
      LINEAR_MIN = 10000
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "flow_accumulation.sgrd",
      FILE = "flow_accumulation.tif"
    )
  )
  
  # Topographic wetness index
  
  slope_r <- paste0(dir, "slope.tif") %>% rast()
  fa_r <- paste0(dir, "flow_accumulation.tif") %>% rast()
  
  s1 <- c(fa_r, slope_r)
  
  get_twi <- function(a, b) {
    b[b == 0] <- 0.001
    out <- log(
      a / (
        tan(b)
      )
    )
    out <- round(out, digits = 2)
    out[!is.finite(out)] <- 0
    return(out)
  }
  
  twi_r <- lapp(s1, fun = get_twi)
  
  writeRaster(
    twi_r,
    datatype = "FLT4S",
    filetype = "GTiff",
    filename = paste0(dir, "TWI.tif"),
    overwrite = TRUE
  )
  
  # Module 15: SAGA wetness index
  
  rsaga.geoprocessor(
    "ta_hydrology",
    module = 15,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      TWI = "TWI_SAGA.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "TWI_SAGA.sgrd",
      FILE = "TWI_SAGA.tif"
    )
  )
  
  # Module 22: LS-Factor
  
  rsaga.geoprocessor(
    "ta_hydrology",
    module = 22,
    env = work_env,
    param = list(
      SLOPE = "slope.sgrd",
      AREA = "flow_accumulation.sgrd",
      LS = "slope_length_factor.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "slope_length_factor.sgrd",
      FILE = "slope_length_factor.tif"
    )
  )
  
  # CHANNELS
  
  # Preprpcess to breach sinks
  
  rsaga.geoprocessor(
    "ta_preprocessor",
    module = 7,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      NOSINKS = "DEM_nosinks.sgrd"
    )
  )
  
  dem_nosinks <- paste0(dir, "DEM_nosinks.sdat") %>% rast()
  
  
  breaches_r <- dem_nosinks < dem_t
  
  # Module 5: Channel Network and Drainage Basins (to create channels)
  
  rsaga.geoprocessor(
    "ta_channels",
    module = 5,
    env = work_env,
    param = list(
      DEM = "DEM_nosinks.sgrd",
      SEGMENTS = "channels_2.shp"
    )
  )
  
  rsaga.geoprocessor(
    "ta_channels",
    module = 5,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      SEGMENTS = "channels.shp"
    )
  )
  
  # Convert channels to grid
  
  rsaga.geoprocessor(
    "grid_gridding",
    module = 0,
    env = work_env,
    param = list(
      INPUT = "channels_2.shp",
      TARGET_TEMPLATE = paste0(DEM, ".sgrd"),
      GRID = "channels_2.sgrd",
      TARGET_DEFINITION = 1
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "channels_2.sgrd",
      FILE = "channels_2.tif"
    )
  )
  
  channels_2_grid <- paste0(dir, "channels_2.tif") %>% rast()
  
  ifel(breaches_r, NA, channels_2_grid * 0 + 1) %>%
    writeRaster(
      filename = paste0(dir, "channels_3.sdat"),
      overwrite = TRUE
    )
  
  rsaga.geoprocessor(
    "grid_gridding",
    module = 0,
    env = work_env,
    param = list(
      INPUT = "channels.shp",
      TARGET_TEMPLATE = paste0(DEM, ".sgrd"),
      GRID = "channels.sgrd",
      TARGET_DEFINITION = 1
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "channels.sgrd",
      FILE = "channels.tif"
    )
  )
  
  # Get local depressions
  
  local_min <- focal(dem_t, fun = "min", na.rm = TRUE) == dem_t
  
  local_min[local_min == 0] <- NA
  
  channels_grid <- paste0(dir, "channels.tif") %>% rast()
  
  channels_plus <- ifel(is.na(channels_grid), local_min, 1)
  
  writeRaster(
    channels_plus,
    filename = paste0(dir, "channels_plus.sdat"),
    overwrite = TRUE
  )
  
  # Horizontal distance to channels
  
  rsaga.geoprocessor(
    "grid_tools",
    module = 10,
    env = work_env,
    param = list(
      SOURCE = "channels.sgrd",
      DISTANCE = "channel_distance_horizontal.sgrd",
      DIST = 500
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "channel_distance_horizontal.sgrd",
      FILE = "channel_distance_horizontal.tif"
    )
  )
  
  # Module 3: Vertical Distance to Channel Network
  
  outnames_ch3 <- c(
    "channel_distance_vertical",
    "channel_base_level"
  )
  
  rsaga.geoprocessor(
    "ta_channels",
    module = 3,
    env = work_env,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      CHANNELS = "channels.sgrd",
      DISTANCE = "channel_distance_vertical.sgrd",
      BASELEVEL = "channel_base_level.sgrd"
    )
  )
  
  outr_ch3 <- outnames_ch3 %>%
    sapply(
      function(x) {
        x %>%
          paste0(dir, ., ".sdat") %>%
          rast()
      }
    ) %>%
    rast()
  
  crs(outr_ch3)
  crs(outr_ch3) <- crs(dem_t)
  setMinMax(outr_ch3)
  minmax(outr_ch3)
  outr_ch3
  
  outr_ch3[[1]][outr_ch3[[1]] < 0] <- 0
  
  outr_ch3 %>%
    math(.,
         "round",
         digits = 2
    ) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = paste0(dir, outnames_ch3, ".tif"),
      overwrite = TRUE
    )
  
  # Downhill gradient to surface water (own calculation)
  
  channel_depth_r <- rast(paste0(dir, "channel_distance_vertical.tif"))
  channel_distance_horizontal_r <- rast(
    paste0(dir, "channel_distance_horizontal.tif")
  )
  
  channel_grad_raw <- atan(channel_depth_r / channel_distance_horizontal_r)
  
  channel_grad <- ifel(channel_distance_horizontal_r == 0, 0, channel_grad_raw)
  
  names(channel_grad) <- "channel_slope_downhill"
  
  channel_grad %>%
    math(.,
         "round",
         digits = 3
    ) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = paste0(dir, "channel_slope_downhill.tif"),
      overwrite = TRUE
    )
  
  # Module 4: Overland Flow Distance to Channel Network
  
  rsaga.geoprocessor(
    "ta_channels",
    module = 4,
    env = work_env,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      CHANNELS = "channels_plus.sgrd",
      DISTANCE = "channel_flow_distance.sgrd",
      DISTVERT = "channel_flow_distance_v.sgrd",
      DISTHORZ = "channel_flow_distance_h.sgrd"
    )
  )
  
  # Downhill flow gradient to surface water (own calculation)
  
  channel_depth_r <- rast(paste0(dir, "channel_flow_distance_v.sdat"))
  channel_distance_horizontal_r <- rast(
    paste0(dir, "channel_flow_distance_h.sdat")
  )
  
  channel_grad_raw <- atan(channel_depth_r / channel_distance_horizontal_r)
  
  channel_grad <- ifel(channel_distance_horizontal_r == 0, 0, channel_grad_raw)
  
  names(channel_grad) <- "channel_flow_slope"
  
  channel_grad %>%
    math(.,
         "round",
         digits = 3
    ) %>%
    writeRaster(
      datatype = "FLT4S",
      filetype = "GTiff",
      filename = paste0(dir, "channel_flow_slope.tif"),
      overwrite = TRUE
    )
  
  # Module 7: Valley depth
  
  rsaga.geoprocessor(
    "ta_channels",
    module = 7,
    env = work_env,
    param = list(
      ELEVATION = paste0(DEM, ".sgrd"),
      VALLEY_DEPTH = "depth_valley_2.sgrd",
      RIDGE_LEVEL = "ridge_level.sgrd"
    )
  )
  
  convertibles <- c("depth_valley_2", "ridge_level")
  
  for (j in 1:length(convertibles)) {
    rsaga.geoprocessor(
      "io_gdal",
      module = 2,
      env = work_env,
      param = list(
        GRIDS = list(paste0(convertibles[j], ".sgrd")),
        FILE = paste0(convertibles[j], ".tif")
      )
    )
  }
  
  # LIGHTING
  
  # Module 5: Topographic Openness
  
  rsaga.geoprocessor(
    "ta_lighting",
    module = 5,
    env = work_env,
    param = list(
      DEM = paste0(DEM, ".sgrd"),
      POS = "openness_positive.sgrd",
      NEG = "openness_negative.sgrd"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "openness_positive.sgrd",
      FILE = "openness_positive.tif"
    )
  )
  
  rsaga.geoprocessor(
    "io_gdal",
    module = 2,
    env = work_env,
    param = list(
      GRIDS = "openness_negative.sgrd",
      FILE = "openness_negative.tif"
    )
  )
  
  # RSAGA TOOLS
  
  mylatitude <- field_polygon %>%
    centroids() %>%
    project("epsg:4326") %>%
    geom() %>%
    as.data.frame() %>%
    .$y
  
  # Potential direct incoming solar radiation
  
  rsaga.pisr2(
    in.dem = paste0(DEM, ".sgrd"),
    out.total.grid = "insolation.sgrd",
    latitude = mylatitude,
    time.step = 1,
    start.date = list(day = 1, month = 12, year = 2020),
    end.date = list(day = 30, month = 11, year = 2021),
    env = work_env
  )
  
  rsaga.geoprocessor("io_gdal",
                     module = 2,
                     env = work_env,
                     param = list(
                       GRIDS = list("insolation.sgrd"),
                       FILE = "insolation.tif"
                     )
  )
  
  # [cleanup] Delete .sgrd files and other SAGA GIS files
  
  patterns <- c("sgrd$", "mgrd$", "prj$", "sdat$", "xml$")
  for (j in 1:length(patterns)) {
    rlist <- list.files(dir, pattern = patterns[j], full.names = TRUE)
    file.remove(rlist)
  }
  
  rm(rlist, j, patterns)
}

# Crop covariates to field extent

library(stringi)
library(tools)

for (i in 1:length(sitenames)) {
  
  dir <- paste0(dir_dat, sitenames[i], "/dem_2m/")
  dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/") %T>%
    dir.create(showWarnings = FALSE)
  
  field_polygon <- paste0(dir_dat, sitenames[i], "/field_polygon.shp") %>%
    vect()
  
  cov_raw <- dir %>%
    list.files(
      pattern = ".tif", full.names = TRUE
    ) %>%
    data.frame(name = .) %>%
    filter(
      !(name %in% grep("channels.tif", name, value = TRUE)),
      !(name %in% grep("channels_2.tif", name, value = TRUE)),
    ) %>%
    unlist() %>%
    unname()
  
  for (j in 1:length(cov_raw)) {
    basename_j <- cov_raw[j] %>% basename() %>% file_path_sans_ext()
    
    outname <- paste0(dir_cov, "/", basename_j, ".tif")
    
    r_masked <- cov_raw[j] %>%
      rast() %>%
      crop(y = field_polygon, mask = TRUE, touches = TRUE)
    
    names(r_masked) <- basename_j
    
    writeRaster(
      r_masked, filename = outname, overwrite = TRUE
    )
  } 
}

i <- 1

dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")

covs_i <- dir_cov %>%
  list.files(full.names = TRUE) %>%
  rast()

plot(covs_i, col = cividis(100), legend = FALSE, axes = FALSE, mar = 0)

# END