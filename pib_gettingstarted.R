# Initial analyses for PIbalance sampling densities

# To do:
# Load data
# Vindum identify and select grid points
# Violin plots for both sites (soil properties)
# Correlation between soil properties and covariates (both sites)
# - preliminary analysis using present data
# - final analysis at 1 m resolution
# New DEM with 1 m accuracy (Vindum + parcels)
# Aggregate rgb to 1 m (both sites)
# Interpolate DUALEM, 1 m (Vindum)
# Covariates at parcel level
# Number og samples vs accuracy
# - Use kmeans
# - only xy
# - using covariates
# - repeat with different folds

library(terra)
library(magrittr)
library(tidyterra)
library(dplyr)
library(ggcorrplot)
library(GGally)
library(PerformanceAnalytics)
library(psych)
library(viridis)

# Load Vindum data
# Soil observations

dir_code <- getwd()
root <- dirname(dir_code)
dir_fig <- paste0(root, "/Figures/") %T>%
  dir.create()

Vindum_obs <- root %>%
  paste0(., "/PIbalance_gis/Vindum/vindum_punkter_top.shp") %>%
  vect() %>%
  filter(
    X %in% seq(0, 1000, by = 20),
    Y %in% seq(0, 1000, by = 20),
  ) %>%
  mutate(
    TOTAL_N = case_when(
      TOTAL_N == "#DIV/0!" ~ NA,
      .default = as.numeric(gsub(",", ".", gsub("\\.", "", TOTAL_N)))
    ),
    ORG_C = case_when(
      ORG_C == "#DIV/0!" ~ NA,
      .default = as.numeric(gsub(",", ".", gsub("\\.", "", ORG_C)))
    )
  ) %>%
  mutate(
    across(
      c(LER, SILT, C_SILT, C_FSAND, FSAND, C_SAND),
      ~ case_when(
        . == 0 ~ NA,
        .default = .
      )
    )
  )

vindum_field <- root %>%
  paste0(., "/PIbalance_gis/Vindum/Vindum_field.shp") %>%
  vect()

Vindum_cov_files <- root %>%
  paste0(., "/Covariates/") %>%
  list.files(full.names = TRUE)

new_ext <- sprc(Vindum_cov_files) %>% ext()

Vindum_cov <- lapply(
  Vindum_cov_files,
  function(x) {
    out <- rast(x)
    out <- extend(out, new_ext)
    return(out)
  }
) %>%
  rast()

plot_ext <- extend(
  ext(Vindum_obs),
  c(0, 250, 0, 0)
)

# Plot observations

tiff(
  paste0(dir_fig, "Vindum_obs_", Sys.Date(), ".tiff"),
  width = 20, height = 10, units = "cm", res = 300
)

Vindum_obs %>%
  select(
    c(HUMUS, LER, SILT, C_SILT, C_FSAND, FSAND, C_SAND, RT, PT, KT, TOTAL_N)
  ) %>%
  plot(
    y = 1:ncol(.),
    legend = "topright",
    mar = c(0, 0, 2.1, 0),
    axes = FALSE,
    box = TRUE,
    ext = plot_ext,
    col = plasma(5)
  )

try(dev.off())
try(dev.off())

# Plot covariates

tiff(
  paste0(dir_fig, "Vindum_cov_", Sys.Date(), ".tiff"),
  width = 20, height = 10, units = "cm", res = 300
)

plot(
  Vindum_cov,
  mar = c(0, 0, 1.5, 4),
  axes = FALSE,
  box = TRUE,
  col = plasma(100),
  nc = 6,
  nr = 4,
  maxnl = 22,
  cex.main = 1
)

try(dev.off())
try(dev.off())

# Extract covariates

vindum_extract <- terra::extract(
  Vindum_cov,
  Vindum_obs,
  ID = FALSE
)

# Correlation plot for soil observations

vindum_obs_sign <- Vindum_obs %>%
  values() %>%
  select(
    c(HUMUS, LER, SILT, C_SILT, C_FSAND, FSAND, C_SAND, RT, PT, KT, TOTAL_N)
  ) %>%
  corr.test()

tiff(
  paste0(dir_fig, "Vindum_cor_obs_", Sys.Date(), ".tiff"),
  width = 16, height = 10, units = "cm", res = 300
)

Vindum_obs %>%
  values() %>%
  select(
    c(HUMUS, LER, SILT, C_SILT, C_FSAND, FSAND, C_SAND, RT, PT, KT, TOTAL_N)
  ) %>%
  cor(use = "pairwise.complete.obs") %>%
  ggcorrplot(
    type = "upper",
    p.mat = vindum_obs_sign$p,
    insig = "blank",
    lab = TRUE,
    lab_size = 3
  )

try(dev.off())
try(dev.off())

# Correlation plot for soil observations vs covariates

tiff(
  paste0(dir_fig, "Vindum_cor_obscov_", Sys.Date(), ".tiff"),
  width = 16, height = 10, units = "cm", res = 300
)

Vindum_obs %>%
  values() %>%
  select(
    c(HUMUS, LER, SILT, C_SILT, C_FSAND, FSAND, C_SAND, RT, PT, KT, TOTAL_N)
  ) %>%
  cor(
    vindum_extract,
    .,
    use = "pairwise.complete.obs"
  ) %>%
  ggcorrplot(method = "square", lab = TRUE, lab_size = 2, tl.cex = 9)

try(dev.off())
try(dev.off())

# END
