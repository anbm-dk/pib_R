# 4: Interpolate gamme ray data

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

# for (i in 1:nrow(sites_df)) {
  # for (i in 3) {
#   
# }

i <- 4

dir_cov <- paste0(dir_dat, sites_df$name[i], "/covariates/")

dem <- dir_cov %>%
  paste0(., "dem_2m.tif") %>%
  rast()
  
gamma <- dir_dat %>%
  paste0(., sitenames[i], "/gamma/concentrations.gpkg") %>%
  vect()

s <- gamma %>%
  select(any_of("Countrate")) %>%
  sf::st_as_sf(.)

variogram <- autofitVariogram(
  as.formula(paste0("Countrate", " ~ 1")),
  s
  # ,
  # fix.values = c(forced_nugget, NA, NA)
)

xy_full <- geom(gamma) %>%
  as.data.frame() %>%
  select(x, y)

gOK <- gstat(
  NULL,
  "Countrate",
  as.formula(paste0("Countrate", " ~ 1")),
  gamma %>%
    select(any_of("Countrate")) %>%
    as.data.frame(geom = "XY"),
  locations = ~ x + y,
  model = variogram$var_model,
  nmin = 6,
  maxdist = 100,
  omax = 2,
  force = TRUE
)

p <- interpolate(dem, gOK, index = 1)
p <- mask(p, dem)

plot(p)


# Boot krig




field_area <- expanse(dem, transform = FALSE) %>%
  select(area) %>%
  unlist() %>%
  unname()

n_pts <- terra::extract(dem, gamma) %>%
  na.omit() %>%
  nrow()

dens_expected <- n_pts / field_area

sigma_pts <- sqrt(field_area / (n_pts * pi))

library(spatstat.geom)
library(spatstat)
dens_out <- ppp(
  xy_full[, 1],
  xy_full[, 2],
  ext(gamma)[1:2],
  ext(gamma)[3:4]
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

for (j in 1:10) {
  set.seed(seeds[j])
  
  bootsample <- gamma %>%
    sample(
      round(sum(weights_dens) / 2),
      prob = weights_dens
    )
  
  id_smalldist <- distance(bootsample) %>%
    apply(1, function(x) {sum(x < (sigma_pts))}) %>%
    unname()
  
  s <- bootsample %>%
    select(any_of("Countrate")) %>%
    sf::st_as_sf(.)
  
  variogram <- autofitVariogram(
    as.formula(paste0("Countrate", " ~ 1")),
    s
    # ,
    # fix.values = c(forced_nugget, NA, NA)
  )
  
  xy <- geom(bootsample) %>%
    as.data.frame() %>%
    select(x, y)
  
  gOK <- gstat(
    NULL,
    "Countrate",
    as.formula(paste0("Countrate", " ~ 1")),
    bootsample %>%
      select(any_of("Countrate")) %>%
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
  
  p <- interpolate(dem, gOK, index = 1)
  p <- mask(p, dem)
  p <- clamp(p, lower = lower_p, upper_p)
  
  plist[[j]] <- p
}
pkrig <- plist %>% rast() %>% mean(na.rm = TRUE)

plot(pkrig)

writeRaster(
  pkrig,
  filename = paste0(dir_cov, "Gamma_Countrate.tiff"),
  overwrite = TRUE,
  names = "Gamma_Countrate"
)

tiff(
  paste0(dir_results, "/Vindum_gamma_map.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

plot(
  pkrig,
  main = "Count rate"
  )

try(dev.off())
try(dev.off())

if (sitenames[i] == "Vindum") {
  obs <- dir_dat %>%
    paste0(., sitenames[i], "/observations.gpkg") %>%
    vect() %>%
    filter(
      Depth == 25,
      is.finite(P)
    )
} else {
  obs <- dir_dat %>%
    paste0(., sitenames[i], "/observations.gpkg") %>%
    vect() %>%
    filter(is.finite(P))
}

obs_extr <- terra::extract(p, obs, bind = TRUE)

head(obs_extr)

cor(obs_extr$Countrate.pred, obs_extr$P)

tiff(
  paste0(dir_results, "/Vindum_gamma_cor.tiff"),
  width = 8, height = 5, units = "cm",
  res = 300
)

ggplot(
  obs_extr,
  aes(x = Countrate.pred, y = P)
) +
  geom_point(shape = 21, bg = "white") +
  ylab("P (mg/g)") +
  xlab("Count rate")

try(dev.off())
try(dev.off())

# Red reflectance for Vindum

red <- paste0(dir_cov, "/ortho_red.tif") %>%
  rast()

tiff(
  paste0(dir_results, "/Vindum_red_ortho_map.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

plot(
  red,
  main = "Red reflectance"
)

try(dev.off())
try(dev.off())

obs_extr2 <- terra::extract(red, obs, bind = TRUE)

head(obs_extr2)

cor(obs_extr2$ortho_red, obs_extr2$P)

tiff(
  paste0(dir_results, "/Vindum_red_cor.tiff"),
  width = 8, height = 5, units = "cm",
  res = 300
)

ggplot(
  obs_extr2,
  aes(x = ortho_red, y = P)
) +
  geom_point(shape = 21, bg = "white") +
  ylab("P (mg/g)") +
  xlab("Red reflectance")

try(dev.off())
try(dev.off())





# xy <- geom(gamma) %>%
#   as.data.frame() %>%
#   select(x, y) %>%
#   as.matrix()
# 
# v <- gamma$Countrate
# 
# tps <- fastTps(xy, v, aRange = 100)
# 
# p <- interpolate(dem, tps)
# p <- mask(p, r)
# plot(p)

# END