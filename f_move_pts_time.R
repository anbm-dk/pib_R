# Function to move DUALEM points depending on measurement times


move_pts_time <- function(
    pts,
    time_shift = 0
) {
  require(MALDIquant)
  
  pts %<>% arrange(GPS_TIME)
  
  pts_xy <- pts %>%
    geom() %>%
    as.data.frame() %>%
    select(x, y) %>%
    as.matrix()
  
  pts_xy_diff <- diff(pts_xy)
  pts_time_diff <- diff(pts$GPS_TIME)
  
  ind_closest <- match.closest(
    pts$GPS_TIME + time_shift,
    pts$GPS_TIME
  )
  
  ind_time_diff <- pts$GPS_TIME + time_shift - pts$GPS_TIME[ind_closest]
  
  neg_time_diff <- as.numeric(ind_time_diff <= 0)
  
  ind_smaller <- ind_closest - neg_time_diff
  
  ind_keep <- ind_smaller > 0
  
  ind_smaller <- ind_smaller[ind_keep]
  
  # print(ind_smaller)
  
  time_residual <- pts$GPS_TIME[ind_keep] + time_shift - pts$GPS_TIME[ind_smaller]
  
  new_pts <- pts_xy %>%
    as.data.frame() %>%
    slice_head(n = nrow(.) - 1) %>%
    rowid_to_column(var = "ID_moved") %>%
    mutate(
      ID_moved = ID_moved + 1,
      x_diff = pts_xy_diff[, 1],
      y_diff = pts_xy_diff[, 2],
      time = pts$GPS_TIME[1:nrow(.)],
      time_diff = pts_time_diff
    )
  
  IDs_orig <- c(1:nrow(pts_xy))[ind_keep]
  
  pts_moved <- new_pts[ind_smaller, ] %>%
    # rowid_to_column(var = "ID_orig") %>%
    mutate(
      ID_orig = IDs_orig,
      time_residual = time_residual,
      x_new = x + x_diff*(time_residual/time_diff),
      y_new = y + y_diff*(time_residual/time_diff),
    ) %>%
    select(-c(x, y)) %>%
    bind_cols(., as.data.frame(pts_xy)[.$ID_orig, ]) %>%
    filter(
      is.finite(x_new),
      is.finite(y_new)
    ) %>%
    rename(
      x_orig = x,
      y_orig = y,
    ) %>%
    mutate(
      xy_gap = sqrt(x_diff^2 + y_diff^2),
      dist_moved = sqrt((x_new - x_orig)^2 + (y_new - y_orig)^2)
    ) %>%
    select(ID_orig, x_orig, y_orig, x_new, y_new, xy_gap, dist_moved)
  
  pts_moved_vect <- pts_moved %>%
    vect(
      geom = c("x_new", "y_new"),
      crs = crs(pts),
      keepgeom = TRUE
    )
  
  values(pts_moved_vect) <- bind_cols(
    values(pts[pts_moved_vect$ID_orig, ]),
    values(pts_moved_vect)
  )
  
  return(pts_moved_vect)
}

# END