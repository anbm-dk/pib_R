# Field-scale validation of DIGIJORD maps

library(terra)
library(devtools)
library(tidyr)
library(tidyterra)
library(magrittr)
library(tibble)
library(dplyr)

library(soilscaler)

my_obs <- list_unwrap(DK_observations, "EPSG:25832")

sitedirs <- c(
  "C:/Users/au542768/Dropbox/AU/P_PIbalance/PIbalance_gis/Field_3",
  "C:/Users/au542768/Dropbox/AU/P_PIbalance/PIbalance_gis/Field_39",
  "C:/Users/au542768/Dropbox/AU/P_PIbalance/PIbalance_gis/Field_47"
)

colnames_keep <- c("Site", "clay", "silt", "sand_f", "sand_c", "SOC")

obs_digijord <- list()

for (i in 1:length(sitedirs)) {
  obs_digijord[[i]] <- paste0(sitedirs[i], "/observations.gpkg") %>%
    vect() %>%
    mutate(
      Site = basename(sitedirs[i]),
      SOC = SOM*0.586,
      clay = Clay*100/(Clay + Silt + Fine_sand + Coarse_sand),
      silt = Silt*100/(Clay + Silt + Fine_sand + Coarse_sand),
      sand_f = Fine_sand*100/(Clay + Silt + Fine_sand + Coarse_sand),
      sand_c = Coarse_sand*100/(Clay + Silt + Fine_sand + Coarse_sand)
    ) %>%
    select(any_of(colnames_keep))
}

obs_digijord

obs_other <- list()

for (i in 1:length(my_obs)) {
  obs_other[[i]] <- my_obs[[i]] %>%
    mutate(
      Site = basename(names(my_obs)[i]),
    ) %>%
    select(any_of(colnames_keep))
}

obs_other

obs_all <- c(obs_digijord, obs_other)

obs_all_merged <- obs_all %>% vect()

# 2014

layers2014 <- c(
  "aclaynor.tif", 
  "asiltnor.tif", 
  "afsandno.tif", 
  "agsandno.tif", 
  "akulstof.tif"
)

rast2014 <- paste0(
  "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/",
  "Texture3D_2014/geotiffs/", layers2014) %>%
  rast()

# 2024

layers2024 <- c(
  "Clay_mean_000_030_cm.tif",
  "Silt_mean_000_030_cm.tif",
  "Fine_sand_mean_000_030_cm.tif",
  "Coarse_sand_mean_000_030_cm.tif"
)

rast2024 <- paste0(
  "O:/Tech_AGRO/Jord/anbm/GIS/Texture_maps_10m/SoilMaps_10m_20240703/",
  "Tekstur2024/depth_000_030_cm/",
  layers2024
) %>%
  c(., paste0(
    "O:/Tech_AGRO/Jord/anbm/GIS/Kulstof2022_usikkerheder_input_parts/",
    "DIGIJORD_SOC_uncertainties/SOC_mean.tif"
  )) %>%
  rast()

extr2014 <- terra::extract(
  rast2014,
  obs_all_merged,
  ID = FALSE
)

extr2024 <- terra::extract(
  rast2024,
  obs_all_merged,
  ID = FALSE
)

results_all <- obs_all_merged %>%
  values() %>%
  rowid_to_column() %>%
  pivot_longer(-c(rowid, Site), values_to = "obs", names_to = "fraction") %>%
  mutate(
    pred2014 = extr2014 %>%
      pivot_longer(
        everything(), 
        values_to = "obs", 
        names_to = "fraction"
        ) %>%
      select(obs) %>%
      unlist(),
    pred2024 = extr2024 %>%
      pivot_longer(
        everything(), 
        values_to = "obs", 
        names_to = "fraction"
      ) %>%
      select(obs) %>%
      unlist()
    ) %>%
  pivot_longer(
    c(pred2014, pred2024),
    values_to = "pred",
    names_to = "map"
  )

summary_all <- results_all %>%
  drop_na() %>%
  group_by(Site, fraction, map) %>%
  summarise(
    R = cor(obs, pred),
    RMSE = sqrt(mean((obs - pred)^2)),
    n = n()
  ) %>%
  filter(
    !(Site %in% names(my_obs)[1:4])
    ) %>%
  arrange(
    map,
    Site,
    match(fraction, colnames_keep)
    ) 

print(summary_all, n = 100)

summary_all %>%
  group_by(fraction, map) %>%
  summarise(R_mean = mean(R),
            RMSE_mean = mean(RMSE)) %>%
  arrange(map, match(fraction, colnames_keep))

obs_all_merged %>%
  geom() %>%
  data.frame() %>%
  mutate(
    Site = obs_all_merged$Site
  ) %>% filter(
    !(Site %in% names(my_obs)[1:4])
  ) %>%
  group_by(Site) %>%
  summarise(
    x = mean(x),
    y = mean(y)
    ) %>%
  vect(geom = c("x", "y"), crs = crs(obs_all_merged), keepgeom = TRUE) %>%
  autoplot()

# END