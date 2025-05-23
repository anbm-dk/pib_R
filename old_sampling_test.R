# Old code


i <- which(sitenames == "Vindum")

dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")

obs <- dir_dat %>%
  paste0(., sitenames[i], "/observations.gpkg") %>%
  vect() %>%
  filter(Depth == 25)

plot(obs, "P", col = cividis(5))

targets <-  c(
  "Clay", "Silt", "Fine_sand", "Coarse_sand",  "SOC",  "Rt", "N", "P", "K"
)

tiff(
  paste0(dir_results, "/Vindum_obs_20240617.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

par(mfrow = c(3,3))

for (i in 1:length(targets)) {
  plot(obs, targets[i], col = cividis(5), main = targets[i], mar = c(1.5, 1.5, 1.5, 3))
}

try(dev.off())
par()

cov_files <- dir_cov %>%
  list.files(pattern = ".tif", full.names = TRUE)

cov_names <- cov_files %>%
  basename() %>%
  file_path_sans_ext() 

cov <- cov_files %>%
  rast()

names(cov) <- cov_names

cov_xy <- cov %>% terra:: subset(c("UTMX", "UTMY"))
cov_selected <- cov %>% terra:: subset(c("ortho_red", "ECa_HCP1m"))

cov_xy_selected <- c(cov_xy, cov_selected)

tiff(
  paste0(dir_results, "/Vindum_cov_examples_20240617.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

plot(cov_xy_selected)

try(dev.off())

obs_cov_all <- terra::extract(cov_selected, obs, bind = TRUE)


# New seed for every iteration

set.seed(1)

obs_cov_val <- obs_cov_all %>%
  slice_sample(prop = 0.25)
obs_cov_train <- obs_cov_all %>%
  filter(!(ID %in% obs_cov_val$ID))

# sample(3:nrow(obs_cov_train), 20)

n_obs_try <- seq.int(
  3,
  min(50, nrow(obs_cov_train)),
  length.out = 20
) %>%
  round(0)

# Matrices for saving results

rmse_xy <- matrix(
  numeric(),
  ncol = length(targets),
  nrow = length(n_obs_try),
  dimnames = list(list(), targets)
)
rmse_cov <- rmse_xy
rmse_xy_cov <- rmse_xy

r2_xy <- rmse_xy
r2_cov <- rmse_xy
r2_xy_cov <- rmse_xy

n_actual_xy <- rmse_xy
n_actual_cov <- rmse_xy
n_actual_xy_cov <- rmse_xy

for (j in 1:length(n_obs_try)) {
  n_obs <- n_obs_try[j]
  
  # k <- 1
  
  # Just xy
  
  set.seed(1)
  
  myclusters_v <- obs_cov_train %>%
    sample_kmeans(input = ., use_xy = TRUE, only_xy = TRUE, clusters = n_obs)
  
  for (k in 1:length(targets)) {
    targ <- targets[k]
    
    traindata <- obs_cov_train[myclusters_v$points$Index, ] %>%
      filter(if_any(targ, is.finite))
    
    if (nrow(traindata) > 2) {
      xydata <- traindata %>%
        select("UTMX", "UTMY") %>%
        values()
      
      targcol <- traindata %>% select(any_of(targ)) %>% values() %>% unlist()
      
      ttt <- train(
        x = xydata,
        y = targcol,
        method = "gaussprPoly",
        trControl = trainControl(method = "LOOCV")
      )
      
      val_p <- values(obs_cov_val) %>%
        select(any_of(names(cov_xy))) %>%
        predict(ttt, .)
      
      val_obs_pred <- obs_cov_val %>%
        select(any_of(targ)) %>%
        values() %>%
        mutate(pred = val_p) %>%
        na.omit()
      
      r2_xy[j, k] <- cor(val_obs_pred[, 1], val_obs_pred[, 2])^2
      rmse_xy[j, k] <- ModelMetrics::rmse(val_obs_pred[, 1], val_obs_pred[, 2])
      n_actual_xy[j, k] <- nrow(traindata)
    }
  }
  
  # With only sensors
  
  set.seed(1)
  
  myclusters_v <- obs_cov_train %>%
    select(names(cov_selected)) %>%
    sample_kmeans(input = ., use_xy = FALSE, clusters = n_obs)
  
  for (k in 1:length(targets)) {
    targ <- targets[k]
    
    traindata <- obs_cov_train[myclusters_v$points$Index, ] %>%
      filter(if_any(targ, is.finite))
    
    if (nrow(traindata) > 2) {
      cov_data <- traindata %>%
        values() %>%
        select(any_of(names(cov_selected)))
      
      targcol <- traindata %>% select(any_of(targ)) %>% values() %>% unlist()
      
      ttt <- train(
        x = cov_data,
        y = targcol,
        method = "gaussprPoly",
        trControl = trainControl(method = "LOOCV")
      )
      
      val_p <- values(obs_cov_val) %>%
        select(any_of(names(cov_selected))) %>%
        predict(ttt, .)
      
      val_obs_pred <- obs_cov_val %>%
        select(any_of(targ)) %>%
        values() %>%
        mutate(pred = val_p) %>%
        na.omit()
      
      r2_cov[j, k] <- cor(val_obs_pred[, 1], val_obs_pred[, 2])^2
      rmse_cov[j, k] <- ModelMetrics::rmse(val_obs_pred[, 1], val_obs_pred[, 2])
      n_actual_cov[j, k] <- nrow(traindata)
    }
  }
  
  # With xy and sensor data
  
  set.seed(1)
  
  myclusters_v <- obs_cov_train %>%
    select(names(cov_selected)) %>%
    sample_kmeans(input = ., use_xy = TRUE, pca = TRUE, clusters = n_obs)
  
  for (k in 1:length(targets)) {
    targ <- targets[k]
    
    traindata <- obs_cov_train[myclusters_v$points$Index, ] %>%
      filter(if_any(targ, is.finite))
    
    if (nrow(traindata) > 2) {
      xydata_plus <- traindata %>%
        values() %>%
        select(any_of(names(cov_xy_selected)))
      
      targcol <- traindata %>% select(any_of(targ)) %>% values() %>% unlist()
      
      ttt <- train(
        x = xydata_plus,
        y = targcol,
        method = "gaussprPoly",
        trControl = trainControl(method = "LOOCV")
      )
      
      val_p <- values(obs_cov_val) %>%
        select(any_of(names(cov_xy_selected))) %>%
        predict(ttt, .)
      
      val_obs_pred <- obs_cov_val %>%
        select(any_of(targ)) %>%
        values() %>%
        mutate(pred = val_p) %>%
        na.omit()
      
      r2_xy_cov[j, k] <- cor(val_obs_pred[, 1], val_obs_pred[, 2])^2
      rmse_xy_cov[j, k] <- ModelMetrics::rmse(val_obs_pred[, 1], val_obs_pred[, 2])
      n_actual_xy_cov[j, k] <- nrow(traindata)
    }
  }
  
  print(paste0("j = ", j))
}
