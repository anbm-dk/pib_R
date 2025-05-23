# Field 3 EMI map


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

set.seed(1234)

seeds <- sample(c(10000:99999), 1000)

dir_dat <- paste0(root, "/PIbalance_gis/")
# sitenames <- c("Field_3", "Field_39", "Field_47", "Trial_070552323", "Vindum")
sitenames <- c("Field_3", "Field_39", "Field_47", "Vindum")

# Do parcels later

sites_df <- data.frame(name = sitenames)

dir_results <- root %>%
  paste0(., "/results/") %T>%
  dir.create()

i <- 1

dir_cov <- paste0(dir_dat, sites_df$name[i], "/covariates/")

dualem_r <- dir_cov %>%
  paste0(., "DUALEM_PRP1m.tiff") %>%
  rast()

tiff(
  paste0(dir_results, "/Field3_EMI_map.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

plot(
  dualem_r,
  main = "ECa (mS/m)"
)

try(dev.off())
try(dev.off())

# END