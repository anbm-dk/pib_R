# Get density weights

get_weightsdens <- function(
    x = NULL,
    rast = NULL,
    dens_mean = NULL,
    sigma = NULL,
    maxw = 1
) {
  require(terra)
  require(magrittr)
  require(spatstat.geom)
  require(spatstat)
  require(dplyr)
  
  if (!is.null(rast)) {
    sum_area <- expanse(rast, transform = FALSE) %>%
      select(area) %>%
      unlist() %>%
      unname()
    
    n_pts <- terra::extract(rast, x) %>%
      na.omit() %>%
      nrow()
    
    ext12 <- ext(rast)[1:2]
    ext34 <- ext(rast)[3:4]
  } else {
    sum_area <- hull(x) %>% expanse()
    
    n_pts <- nrow(x)
    
    ext12 <- ext(x)[1:2]
    ext34 <- ext(x)[3:4]
  }
  
  if (is.null(dens_mean)) { dens_mean <- n_pts / sum_area }
  if (is.null(sigma)) { sigma <- sqrt(sum_area / (n_pts * pi)) }
  
  dens_out <- ppp(
    geom(x) %>%
      as.data.frame() %>%
      select("x") %>% 
      unlist() %>% 
      unname(),
    geom(x) %>%
      as.data.frame() %>%
      select("y") %>% 
      unlist() %>% 
      unname(),
    ext12,
    ext34
  ) %>%
    density(
      sigma = sigma,
      at = "points",
      leaveoneout = FALSE
    )
  
  attributes(dens_out) <- NULL
  
  weights_dens <- dens_mean / dens_out
  
  if (!is.null(maxw)) {
    weights_dens[weights_dens > maxw] <- maxw
  }
  
  print(dens_mean)
  print(sigma)
  print(sum_area)
  
  return(weights_dens)
}

# library(terra)
# library(magrittr)
# library(sp)
# library(sf)
# library(obliquer)
# 
# v <- unwrap(Vindum_SOM)
# myrast <- unwrap(Vindum_covariates)[[1]]
# 
# myweights <- get_weightsdens(v, rast = myrast)
# 
# v2 <- v
# v2$w <- myweights
# 
# plot(v2, "w")
