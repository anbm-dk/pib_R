# 5: Analyses for PIbalance

# TO DO:
# Weighted sampling for k-means.
# Wighted accuracy assessment.
# Move to stationary PC for efficiency.

# library(devtools)
# install_github("anbm-dk/samplekmeans", force = TRUE)
library(samplekmeans)
library(terra)
library(magrittr)
library(viridis)
library(tidyr)
library(tidyterra)
library(tools)
library(fields) #tps
library(dplyr)
library(caret)
library(tibble)
library(corrplot)

if (!exists("started")) {
  wd <- getwd()
  setwd("..")
  root <- getwd()
  setwd(wd)
  rm(wd)
  started <- TRUE
}

dir_dat <- paste0(root, "/PIbalance_gis/")
# sitenames <- c("Field_3", "Field_39", "Field_47", "Trial_070552323", "Vindum")

sitenames <- c("Field_3", "Field_39", "Field_47", "Vindum")

dir_results <- root %>%
  paste0(., "/results/") %T>%
  dir.create()


# Create summary and figures

tiff(
  paste0(dir_results, "/fields_overview_v1.tiff"),
  width = 24, height = 16, units = "cm",
  res = 300
)

par(mfrow = c(2, 2))

for(i in 1:length(sitenames)) {
  dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")
  
  if (sitenames[i] == "Vindum") {
    obs <- dir_dat %>%
      paste0(., sitenames[i], "/observations.gpkg") %>%
      vect() %>%
      filter(Depth == 25)
  } else {
    obs <- dir_dat %>%
      paste0(., sitenames[i], "/observations.gpkg") %>%
      vect()
  }
  
  field_boundary <- dir_dat %>%
    paste0(., sitenames[i], "/field_polygon.shp") %>%
    vect()
  
  dem <- dir_cov %>%
    paste0(., "dem_2m.tif") %>%
    rast()
  
  plot(
    dem,
    main = sitenames[i],
    plg = list(title="Elevation (m)")
    )
  
  plot(field_boundary, add = TRUE)
  
  # plot(obs, "P", col = cividis(5))
}

try(dev.off())
try(dev.off())

par()

# Summary statistics

areas <- numeric()

obs_all <- list()

for(i in 1:length(sitenames)) {
  dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")
  
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
  
  field_boundary <- dir_dat %>%
    paste0(., sitenames[i], "/field_polygon.shp") %>%
    vect()
  
  areas[i] <- terra::expanse(field_boundary, unit = "ha")
  
  obs_values <- obs %>%
    values %>%
    mutate(
    Site = sitenames[i]
  )
  
  obs_all[[i]] <- obs_values
}

areas

summaries_sites <- obs_all %>%
  bind_rows() %>%
  group_by(Site) %>%
  summarise(across(
    .cols = where(is.numeric), 
    .fns = list(Mean = mean, SD = sd), na.rm = TRUE, 
    .names = "{col}_{fn}"
  ))

summaries_sites$area <- areas

sites_n <- obs_all %>%
  bind_rows() %>%
  group_by(Site) %>%
  summarise(n = n())

sites_n

summaries_sites$n <- sites_n$n

write.table(
  summaries_sites,
  paste0(dir_results, "/site_summary_v1.csv"),
  sep = ";",
  row.names = FALSE
  )


# Correlation with P

cor_list <- list()
covobs_list <- list()

for(i in 1:length(sitenames)) {
  dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")
  
  covariates <- dir_cov %>% list.files(".tif", full.names = TRUE) %>% rast()
  
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
  
  obs_cov <- terra::extract(covariates, obs, bind = TRUE) %>%
    values() %>%
    select(any_of(c("P", names(covariates)))) %>%
    drop_na() %>%
    mutate(Site = sitenames[i])
  
  cor_list[[i]] <- cor(obs_cov$P, obs_cov[, c(-1, -ncol(obs_cov))]) %>%
    as.data.frame()
  
  covobs_list[[i]] <- obs_cov
}

cor_list2 <- cor_list %>% bind_rows() %>% t()

colnames(cor_list2) <- sitenames

cor_list2 %<>%
  as.data.frame() %>%
  rownames_to_column() %>%
  drop_na() %>%
  rowwise() %>%
  mutate(
    mean = mean(c(Field_3, Field_39, Field_47, Vindum)),
    mean_2 = mean(abs(c(Field_3, Field_39, Field_47, Vindum))),
    abs_mean = abs(mean)
  )

cor_list2 %>%
  filter(rowname == "mean")

cor_list2 %>%
  arrange(-abs_mean) %>%
  as.data.frame()

cor_list2 %>%
  arrange(-mean_2) %>%
  as.data.frame()

tiff(
  paste0(dir_results, "/ECa_cor_all4.tiff"),
  width = 12, height = 7.5, units = "cm",
  res = 300
)

covobs_list %>%
  bind_rows() %>%
  ggplot(
    aes(x = mean, y = P)
  ) +
  geom_point(shape = 21, bg = "white") +
  facet_wrap(~ Site, scales = "free", nrow = 2) +
  xlab("ECa (mS/m)") +
  ylab("P (mg/g)")

try(dev.off())
try(dev.off())


vars_DSM10m <- c("clay_10m", "silt_10m", "coarse_sand_10m", "SOC_10m")

vars_TOPO_top10 <- c(
  "position_topographic_multi", "convergence_radius", "sinks", "height_slope", 
  "MRVBF", "depth_valley_2", "height_normalized", "height_standardized", "TWI",
  "channel_flow_slope"
  )

covcor_list <- covobs_list %>%
  lapply(
    function(x) {
      out <- x %>%
        select(any_of(cor_list2$rowname)) %>%
        select(any_of(vars_TOPO_top10)) %>%
        select(-any_of(c("curvature_flow_line", "curvature_total"))) %>%
        cor()
      return(out)
    }
  )

covcor_mean <- Reduce("+", covcor_list) / length(covcor_list)

covcor_mean %>%
  corrplot(type = "upper")

covcor_mean %>%
  as.data.frame() %>%
  arrange(-abs(position_topographic_multi))

covcor_mean %>%
  as.data.frame() %>%
  mutate(mean = abs(MRVBF) + abs(position_topographic_multi)) %>%
  arrange(-abs(mean))

covcor_mean %>%
  as.data.frame() %>%
  mutate(mean = abs(MRVBF) + abs(position_topographic_multi) + abs(sinks)) %>%
  arrange(-abs(mean))

vars_TOPO <- c(
  "position_topographic_multi", "MRVBF", "sinks", "TWI"
)

covariates %>%
  select(any_of(vars_DSM10m)) %>%
  plot()

covariates %>%
  select(any_of(vars_TOPO)) %>%
  plot()

vars_xy <- c("UTMX", "UTMY")

varlist <- list(
  vars_xy, vars_DSM10m, vars_TOPO
)

names_cat <- c("XY", "DSM10m", "TOPO")

var_indices <- expand.grid(
  c(NA, 1),
  c(NA, 2),
  c(NA, 3)
) %>%
  apply(
    ., 1, 
        function(x) {
          out <- unname(x[is.finite(x)])
          return(out)
        }
    )

var_indices <- var_indices[lengths(var_indices) != 0]

var_indices

names_cat_list <- lapply(
  var_indices,
  function (x) {
    out <- names_cat[x] %>%
      paste(collapse = "+")
    return(out)
  }
  ) %>% unlist()

names_cat_list

var_combination_list <- lapply(
  var_indices,
  function (x) {
    out <- varlist %>%
      magrittr::extract(x) %>%
      unlist()
    return(out)
  }
)

# Priorities
# a: Focus on P.
# b: Realistic sampling densities.
# c: Correlation with other soil properties (4 big fields).
# d: Correlation with topography (4 big fields).
# e: Variograms (all fields)
# f: Sensors:
#     - EMI
#     - RGB
#     - Gamma
# g: Topographic variables.
# h: 10 m maps.


# Later:
# - Revise k-means sampling function (candidate points, function for dataframes, feature weighting)


# Still missing:
# 1: Crop topographic covariates to field extent. [ok]
# 2: Resample 10 texture maps to field level [ok]
# 3: Krige the DUALEM measurements from Vindum to 2 m resolution [ok]
# 4: Do the same for the other sites when the measurements are ready.
# 5: Add xy rasters to the covariates [ok]
# 6: Try different numbers of samples:
## 6.1: Split training/test data. Repeated splits to ensure robustness.
# Sampling
## 6.3: K-means samples
# Regression:
## 6.4: IDW interpolation

# Discarded/suspended ideas:
## 6.2: Random samples
## 6.5: Cluster-based regression?
## 6.6: Linear regression
## 6.7: Baseline: Mean of samples (start with just one sample).

# First test, just for Vindum:
# 1: Covariates: A: xy only, B: EC + red band, C: A + B
# 2: Sampling: k-means
# 3: Regression: thin plate spline
# Use same points for all targets

time0 <- Sys.time()

n_min <- 5
sample_dens_max <- 2

list_results <- list()

set.seed(1)

seeds <- sample(c(1000:10000), 50)
n_sample_sizes_rep <- 10

# for (i in 1:length(sitenames)) {
#   
#   dir_cov <- paste0(dir_dat, sitenames[i], "/covariates/")
#   
#   covariates <- dir_cov %>% list.files(".tif", full.names = TRUE) %>% rast()
#   
#   if (sitenames[i] == "Vindum") {
#     obs <- dir_dat %>%
#       paste0(., sitenames[i], "/observations.gpkg") %>%
#       vect() %>%
#       filter(
#         Depth == 25,
#         is.finite(P)
#       )
#   } else {
#     obs <- dir_dat %>%
#       paste0(., sitenames[i], "/observations.gpkg") %>%
#       vect() %>%
#       filter(is.finite(P))
#   }
#   
#   obs %<>% select(-any_of(c("UTMX", "UTMY")))
#   
#   cov_files <- dir_cov %>%
#     list.files(pattern = ".tif", full.names = TRUE)
#   
#   cov_names <- cov_files %>%
#     basename() %>%
#     file_path_sans_ext() 
#   
#   names(covariates) <- cov_names
#   
#   obs_cov <- terra::extract(covariates, obs, bind = TRUE) %>%
#     select(any_of(c("ID", "P", names(covariates)))) %>%
#     drop_na()
#   
#   field_boundary <- dir_dat %>%
#     paste0(., sitenames[i], "/field_polygon.shp") %>%
#     vect()
#   
#   area_ha <- terra::expanse(field_boundary, unit = "ha")
#   
#   val_prop <- max(c(0.25, 10/nrow(obs_cov)))
#   
#   n_max <- min(
#     c(floor(nrow(obs_cov)*(0.5)),
#       area_ha*sample_dens_max
#     )
#   )
#   
#   for(j in 1:length(seeds)) {
#     
#     set.seed(seeds[j])
#     
#     ind_val <- sample(1:nrow(obs_cov), val_prop*nrow(obs_cov))
#     
#     obs_cov$ID_new_v2 <- 1:nrow(obs_cov)
#     
#     obs_cov_val <- obs_cov[ind_val, ]
#     obs_cov_train <- obs_cov[-ind_val, ]
#     
#     ns_j <- sample(c(n_min:n_max), n_sample_sizes_rep) %>%
#       sort()
#     
#     for (s in 1:length(ns_j)) {
#       n_j <- ns_j[s]
#       
#       for (k in 1:length(names_cat_list)) {
#         
#         if (names_cat_list[k] == "XY") {
#           use_xy_only <- TRUE
#           use_coordinates <- TRUE
#         } else {
#           use_xy_only <- FALSE
#           use_coordinates <- sum(vars_xy %in% var_combination_list[[k]]) > 0
#         }
#         
#         set.seed(seeds[j])
#         
#         myclusters_v <- obs_cov_train %>%
#           select(any_of(var_combination_list[[k]])) %>%
#           sample_kmeans(
#             input = .,
#             clusters = n_j,
#             only_xy = use_xy_only,
#             use_xy = use_coordinates,
#             pca = !use_xy_only
#           )
#         
#         traindata <- obs_cov_train[myclusters_v$points$Index, ] %>%
#           values() %>%
#           select(any_of(c("P", var_combination_list[[k]])))
#         
#         upperbound <- mean(traindata$P) + sd(traindata$P)*3
#         
#         if (use_xy_only) {
#           set.seed(seeds[j])
#           
#           vars_use <- vars_xy
#           
#           ttt <- train(
#             as.formula(
#               paste("P ~",
#                     paste(
#                       vars_use,
#                       collapse = "+"
#                       )
#                     )
#               ),
#             traindata,
#             method = "gaussprPoly",
#             trControl = trainControl(
#               method = "cv",
#               number = 5,
#               predictionBounds = c(0, upperbound),
#             )
#           )
#         } else {
#           set.seed(seeds[j])
#           
#           ttt0 <- train(
#             as.formula(
#               paste("P ~",
#                     paste(
#                       var_combination_list[[k]],
#                       collapse = "+")
#                     )
#               ),
#             traindata,
#             method = "gaussprPoly",
#             preProcess = c(
#               "nzv"
#               # ,
#               # "pca"
#               ),
#             trControl = trainControl(
#               method = "cv",
#               number = 5,
#               predictionBounds = c(0, upperbound),
#               preProcOptions = list(
#                 freqCut = 80/20,
#                 uniqueCut = 10
#               )
#             )
#           )
#           
#           max_covariates <- min(
#               nrow(traindata) - 2,
#               length(var_combination_list[[k]])
#             )
#           
#           if (use_coordinates) {
#             min_covariates <- 3
#             
#             vars_ordered <- varImp(ttt0)$importance %>%
#               arrange(-Overall) %>%
#               rownames_to_column() %>%
#               filter(!(rowname %in% vars_xy)) %>%
#               select(rowname) %>%
#               unlist() %>%
#               unname() %>%
#               c(vars_xy, .)
#             
#           } else {
#             min_covariates <- 1
#             
#             vars_ordered <- varImp(ttt0)$importance %>%
#               arrange(-Overall) %>%
#               rownames_to_column() %>%
#               select(rowname) %>%
#               unlist() %>%
#               unname()
#           }
#           
#           
#           nvars <- c(min_covariates:max_covariates)
#           
#           rmses_vars <- sapply(
#             nvars,
#             function (x) {
#               vars_use <- vars_ordered[1:x]
#               
#               if (length(vars_use) == 1) {
#                 preprocess_vars <- NULL
#               } else {
#                 preprocess_vars <- c(
#                   "nzv"
#                   # ,
#                   # "pca"
#                   )
#               }
#               
#               set.seed(seeds[j])
#               
#               tttx <- train(
#                 as.formula(paste("P ~",
#                                  paste(vars_use, collapse = "+"))),
#                 traindata,
#                 method = "gaussprPoly",
#                 preProcess = preprocess_vars,
#                 trControl = trainControl(
#                   method = "cv",
#                   number = 5,
#                   predictionBounds = c(0, upperbound),
#                   preProcOptions = list(
#                     freqCut = 80/20,
#                     uniqueCut = 10
#                   )
#                 )
#               )
#               
#               out <- tttx$results$RMSE %>%
#                 min(na.rm = TRUE)
#               
#               return(out)
#             }
#           )
#           
#           vars_use <- vars_ordered[1:nvars[which.min(rmses_vars)]]
#           
#           if (length(vars_use) == 1) {
#             preprocess_vars <- NULL
#           } else {
#             preprocess_vars <- c("nzv", "pca")
#           }
#           
#           set.seed(seeds[j])
#           
#           ttt <- train(
#             as.formula(paste("P ~",
#                              paste(vars_use, collapse = "+"))),
#             traindata,
#             method = "gaussprPoly",
#             preProcess = preprocess_vars,
#             trControl = trainControl(
#               method = "cv",
#               number = 5,
#               predictionBounds = c(0, upperbound),
#               preProcOptions = list(
#                 freqCut = 80/20,
#                 uniqueCut = 10
#               )
#             )
#           )
#           
#         }
#         
#         
#         val_p <- values(obs_cov_val) %>%
#           predict(ttt, .)
#         
#         val_obs_pred <- obs_cov_val %>%
#           select(P) %>%
#           values() %>%
#           mutate(pred = val_p) %>%
#           na.omit()
#         
#         n_per_ha <- nrow(traindata)/area_ha
#         
#         results_k <- data.frame(
#           Site = sitenames[i],
#           Rep = j,
#           Seed = seeds[j],
#           n = nrow(traindata),
#           n_per_ha = n_per_ha,
#           R2 = cor(val_obs_pred[, 1], val_obs_pred[, 2])^2,
#           RMSE = ModelMetrics::rmse(val_obs_pred[, 1], val_obs_pred[, 2]),
#           Layers = names_cat_list[k],
#           vars_used = paste0(vars_use, collapse = "+")
#         )
#         
#         list_results[[length(list_results) + 1]] <- results_k
#         
#         print(results_k)
#       }
#     }
#   }
# }
# 
# all_results <- bind_rows(list_results)

# write.table(
#   all_results,
#   paste0(dir_results, "/sampletest_fourfields_P.csv") ,
#   sep = ";",
#   row.names = FALSE
# )

# time1 <- Sys.time()
# 
# time1 - time0

all_results <- read.csv(
  paste0(dir_results, "/sampletest_fourfields_P.csv"),
    sep = ";"
)

# Sampling densities

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  group_by(Site) %>%
  summarise(
    mindens = min(n_per_ha),
    maxdens = max(n_per_ha)
  )

mindens_t <- 0.43
maxdens_t <- 0.72

# RMSE per method, plot for each site, line only (effect of n samples is small)

tiff(
  paste0(dir_results, "/RMSE_Site.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  ggplot(aes(x = n_per_ha, y = RMSE, col = Layers)) +
  facet_wrap(~ Site, scales = "free", nrow = 2) +
  geom_point() +
  geom_smooth(
    # method = "lm",
    # span = 1,
    se = FALSE
    ) +
  xlab("n/ha")
# +
#   ylim(0, NA)

try(dev.off())
try(dev.off())

# R2 per method, plot for each site, line only

tiff(
  paste0(dir_results, "/R2_Site.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)


all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  ggplot(aes(x = n_per_ha, y = R2, col = Layers)) +
  facet_wrap(~ Site, scales = "free", nrow = 2) +
  geom_point() +
  geom_smooth(
    # method = "lm",
    # span = 1,
    se = FALSE
  ) +
  ylim(0, NA) +
  ylab(bquote(R^2)) +
  xlab("n/ha")

try(dev.off())
try(dev.off())

# Violin plot for RMSE, 0.43 - 0.72 samples per ha

tiff(
  paste0(dir_results, "/Violin_RMSE.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  # filter(RMSE < 2) %>%
  filter(
    n_per_ha < maxdens_t,
    n_per_ha > mindens_t
    ) %>%
  ggplot(aes(x = Layers, y = RMSE, fill = Layers)) +
  geom_violin() +
  facet_wrap(~ Site, nrow = 2, scales = "free_x") +
  stat_summary(fun.y = mean, geom="point", size = 2) +
  coord_flip() +
  xlab(NULL)
# +
#   ylim(0, NA)

try(dev.off())
try(dev.off())

# Violin plot for R2, 0.43 - 0.72 samples per ha

tiff(
  paste0(dir_results, "/Violin_R2.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  filter(RMSE < 2) %>%
  filter(
    n_per_ha < maxdens_t,
    n_per_ha > mindens_t
  ) %>%
  ggplot(aes(x = Layers, y = R2, fill = Layers)) +
  geom_violin(, scale = "width", bounds = c(0, Inf)) +
  facet_wrap(~ Site, nrow = 2, scales = "free_x") +
  stat_summary(fun.y = mean, geom="point", size = 2) +
  coord_flip() +
  ylab(bquote(R^2)) +
  xlab(NULL)

try(dev.off())
try(dev.off())


# +
#   ylim(0, NA)

# Calculate standardized accuracies, by field

acc_means_Sites <- all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  filter(
    n_per_ha < maxdens_t,
    n_per_ha > mindens_t
  ) %>%
  group_by(Site) %>%
  reframe(
    R2_Site = mean(R2),
    RMSE_Site = mean(RMSE)
  )

acc_means_Sites

acc_means_all <- acc_means_Sites %>%
  reframe(
    R2_all = mean(R2_Site),
    RMSE_all = mean(RMSE_Site)
  ) %>%
  unlist()

acc_means_all

# Standardized accuracies, by method

acc_means_Layers <- all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  filter(
    n_per_ha < maxdens_t,
    n_per_ha > mindens_t
  ) %>%
  group_by(Site, Layers) %>%
  reframe(
    R2_mean = mean(R2),
    RMSE_mean = mean(RMSE)
  ) %>%
  group_by(Layers) %>%
  reframe(
    R2_lyr = mean(R2_mean),
    RMSE_lyr = mean(RMSE_mean)
  )

acc_means_Layers

acc_means_Layers_Sites <- all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  filter(
    n_per_ha < maxdens_t,
    n_per_ha > mindens_t
  ) %>%
  group_by(Site, Layers) %>%
  reframe(
    R2_lyrSite = mean(R2),
    RMSE_lyrSite = mean(RMSE)
  )

acc_means_Layers_Sites

acc_summary_scalers <- acc_means_Layers_Sites %>%
  left_join(., acc_means_Sites, by = "Site") %>%
  left_join(., acc_means_Layers, by = "Layers") %>%
  mutate(
    # R2_Site_mult = R2_Site / acc_means_all[1],
    # RMSE_Site_mult = RMSE_Site / acc_means_all[2],
    R2_mult = R2_lyr / R2_lyrSite,
    RMSE_mult = RMSE_lyr / RMSE_lyrSite
    # R2_mult = R2_Site_mult * R2_lyr_mult,
    # RMSE_mult = RMSE_Site_mult * RMSE_lyr_mult
  )

acc_summary_scalers %>% as.data.frame()

# Mean for each field should be 1
# Mean for each method should be the mean ACC for the method divided by the
# overall mean acc

# Standardized R2, one plot

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
      )
  ) %>%
  right_join(., acc_means, by = "Site") %>%
  mutate(
    R2_scaled = R2 / R2_mean,
    RMSE_scaled = RMSE / RMSE_mean
  ) %>%
  ggplot(aes(x = n_per_ha, y = R2_scaled, col = Layers)) +
  # geom_point() +
  geom_smooth(
    # method = "lm"
    # se = FALSE
  ) +
  ylim(0, NA) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed")

# Standardized RMSE, one plot (kind of flat)

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_means, by = "Site") %>%
  mutate(
    R2_scaled = R2 / R2_mean,
    RMSE_scaled = RMSE / RMSE_mean
  ) %>%
  ggplot(aes(x = n_per_ha, y = RMSE_scaled, col = Layers)) +
  # geom_point() +
  geom_smooth(
    # method = "lm"
    # se = FALSE
  ) +
  # ylim(0, NA) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed")

# Standardized RMSE, by method, with points

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_means, by = "Site") %>%
  mutate(
    R2_scaled = R2 / R2_mean,
    RMSE_scaled = RMSE / RMSE_mean
  ) %>%
  ggplot(aes(x = n_per_ha, y = RMSE_scaled, col = Layers)) +
  facet_wrap(~ Layers, nrow = 2) +
  geom_point() +
  geom_smooth(
    col = "black"
    # method = "lm"
    # se = FALSE
  ) +
  # ylim(0, NA) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed")

# Standardized R2, by method, with points

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_means, by = "Site") %>%
  mutate(
    R2_scaled = R2 / R2_mean,
    RMSE_scaled = RMSE / RMSE_mean
  ) %>%
  ggplot(aes(x = n_per_ha, y = R2_scaled, col = Layers)) +
  facet_wrap(~ Layers, nrow = 2) +
  geom_point() +
  geom_smooth(
    col = "black"
    # method = "lm"
    # se = FALSE
  ) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed")

# Using summary scalers

# Standardized RMSE, using summary scalers, one plot

tiff(
  paste0(dir_results, "/RMSE_standardized.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_summary_scalers, by = c("Site", "Layers")) %>%
  mutate(
    R2_scaled = R2*R2_mult,
    RMSE_scaled = RMSE*RMSE_mult
  ) %>%
  ggplot(aes(x = n_per_ha, y = RMSE_scaled, col = Layers)) +
  # geom_point() +
  geom_smooth(
    # formula= (y ~ log(x)),
    # method = "lm",
    se = FALSE
  ) +
  # ylim(0, NA) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed") +
  geom_vline(xintercept = 1) +
  xlab("n/ha") +
  ylab("Scaled RMSE")

try(dev.off())
try(dev.off())

# Standardized R2, using summary scalers, one plot

tiff(
  paste0(dir_results, "/R2_standardized.tiff"),
  width = 16, height = 10, units = "cm",
  res = 300
)

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_summary_scalers, by = c("Site", "Layers")) %>%
  mutate(
    R2_scaled = R2*R2_mult,
    RMSE_scaled = RMSE*RMSE_mult
  ) %>%
  ggplot(aes(x = n_per_ha, y = R2_scaled, col = Layers)) +
  # geom_point() +
  geom_smooth(
    # formula = (y ~ log(x)),
    # method = "lm"
    # se = FALSE
  ) +
  # ylim(0, NA) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed") +
  geom_vline(xintercept = 1) +
  ylab(bquote("Scaled"~R^2)) +
  xlab("n/ha")

try(dev.off())
try(dev.off())

# Standardized RMSE, by method, with points, using summary scalers

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_summary_scalers, by = c("Site", "Layers")) %>%
  mutate(
    R2_scaled = R2*R2_mult,
    RMSE_scaled = RMSE*RMSE_mult
  ) %>%
  ggplot(aes(x = log(n_per_ha), y = log(RMSE_scaled), col = Layers)) +
  facet_wrap(~ Layers, nrow = 2) +
  geom_point() +
  geom_smooth(
    col = "black",
    # formula= (y ~ log(x)),
    # formula = y ~ log(x),
    method = "lm"
    # se = FALSE
  ) +
  ylim(0, NA) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed") +
  xlab("n/ha")

# Standardized R2, by method, with points, using summary scalers

all_results %>%
  mutate(
    R2 = case_when(
      is.na(R2) ~ 0,
      .default = R2
    )
  ) %>%
  right_join(., acc_summary_scalers, by = c("Site", "Layers")) %>%
  mutate(
    R2_scaled = R2*R2_mult,
    RMSE_scaled = RMSE*RMSE_mult
  ) %>%
  ggplot(aes(x = n_per_ha, y = R2_scaled, col = Layers)) +
  facet_wrap(~ Layers, nrow = 2) +
  geom_point() +
  geom_smooth(
    col = "black"
    # formula= (y ~ log(x)),
    # method = "lm"
    # se = FALSE
  ) +
  geom_vline(xintercept = mindens_t, linetype = "dashed") +
  geom_vline(xintercept = maxdens_t, linetype = "dashed") +
  geom_vline(xintercept = 1) +
  ylab(bquote(R^2)) +
  xlab("n/ha")

# To do:
# Make maps for XY+TOPO, 1 sample per 2 ha, for each site.
# Find the most important variables + effect of sample density

# Map from last model

pred_r <- predict(covariates, ttt, na.rm = TRUE)

plot(
  pred_r,
  main = "last model"
)

# Example maps

my_method <- "gam"

ttt2 <- train(
  as.formula(paste("P ~",
                   paste(vars_use, collapse = "+"))),
  traindata,
  method = my_method,
  preProcess = "pca",
  trControl = trainControl(method = "cv", number = 5)
)

ttt2

pred_r <- predict(covariates, ttt, na.rm = TRUE)

plot(
  pred_r,
  main = my_method
  )



# New features for k-means sampler:
## 1: Candidate points (clusters based on raster, only select points from set)
## 2: Output distance layers?
## 3: Cluster-based regression?
## 4: Sub-clusters? (n samples per cluster)

# Figure this out:
## How to select variables (clustering, regression)
### - Correlation from other fields?
### - Feature selection?
### - A priori selection?
### - Variable number based on the number of observations?
## Best setting for individual targets and/or all targets?
## Multivariate regression?

# https://www.rdocumentation.org/packages/fields/versions/15.2/topics/Tps 

# Multivariate Adaptive Regression Splines	gcvEarth
# lmStepAIC
# gaussprRadial
# Projection Pursuit Regression
# Independent Component Regression	icr
# Partial Least Squares	kernelpls
# Gaussian Process gaussprLinear
# Principal Component Analysis	pcr
# Multivariate Adaptive Regression Splines	gcvEarth
# Partial Least Squares	pls

# END