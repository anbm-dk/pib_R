# Weighted summary functions

# Weighted RMSE
get_RMSEw <- function(d, w) {
  if (nrow(d) == 0) {
    out <- NA
  } else {
    sqe <- w * (d[, 1] - d[, 2])^2
    msqe <- sum(sqe) / sum(w)
    out <- sqrt(msqe)
  }
  return(out)
}

# Weighted R^2
get_R2w <- function(d, w) {
  require(boot)
  if (nrow(d) < 3) {
    out <- NA
  } else {
    require(boot)
    out <- boot::corr(d[, 1:2], w)^2
  }
  return(out)
}

# Weighted summary function
WeightedSummary <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  get_RMSEw <- function(d, w) {
    if (nrow(d) == 0) {
      out <- NA
    } else {
      sqe <- w * (d[, 1] - d[, 2])^2
      msqe <- sum(sqe) / sum(w)
      out <- sqrt(msqe)
    }
    return(out)
  }
  get_R2w <- function(d, w) {
    require(boot)
    if (nrow(d) < 3) {
      out <- NA
    } else {
      require(boot)
      out <- boot::corr(d[, 1:2], w)^2
    }
    return(out)
  }
  out <- numeric()
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c("RMSEw", "R2w")
  return(out)
}

# Weighted summary function with log transformation
WeightedSummary_log <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  get_RMSEw <- function(d, w) {
    if (nrow(d) == 0) {
      out <- NA
    } else {
      sqe <- w * (d[, 1] - d[, 2])^2
      msqe <- sum(sqe) / sum(w)
      out <- sqrt(msqe)
    }
    return(out)
  }
  get_R2w <- function(d, w) {
    require(boot)
    if (nrow(d) < 3) {
      out <- NA
    } else {
      require(boot)
      out <- boot::corr(d[, 1:2], w)^2
    }
    return(out)
  }
  out <- numeric()
  data[, 1:2] <- log(data[, 1:2])
  data <- data[is.finite(rowSums(data)), ]
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c("RMSEw_log", "R2w_log")
  return(out)
}

# Weighted summary function with square root transformation
WeightedSummary_sqrt <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  get_RMSEw <- function(d, w) {
    if (nrow(d) == 0) {
      out <- NA
    } else {
      sqe <- w * (d[, 1] - d[, 2])^2
      msqe <- sum(sqe) / sum(w)
      out <- sqrt(msqe)
    }
    return(out)
  }
  get_R2w <- function(d, w) {
    require(boot)
    if (nrow(d) < 3) {
      out <- NA
    } else {
      out <- boot::corr(d[, 1:2], w)^2
    }
    return(out)
  }
  out <- numeric()
  data[, 1:2] <- sqrt(data[, 1:2])
  data <- data[is.finite(rowSums(data)), ]
  out[1] <- get_RMSEw(data[, 1:2], data$weights)
  out[2] <- get_R2w(data[, 1:2], data$weights)
  names(out) <- c("RMSEw_sqrt", "R2w_sqrt")
  return(out)
}

# Weighted summary for drainage classes as numeric vector

WeightedSummary_DCnum <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  get_MAEw <- function(d, w) {
    require(stats)
    if (nrow(d) == 0) {
      out <- NA
    } else {
      ae <- abs(d[, 1] - d[, 2])
      out <- stats::weighted.mean(x = ae, w = w, na.rm = TRUE)
    }
    return(out)
  }
  get_OAw <- function(d, w) {
    require(stats)
    if (nrow(d) == 0) {
      out <- NA
    } else {
      d_int <- round(d, digits = 0)
      d_equal <- d_int[, 1] == d_int[, 2]
      out <- stats::weighted.mean(x = d_equal, w = w, na.rm = TRUE)
    }
    return(out)
  }
  out <- numeric()
  out[1] <- get_MAEw(data[, 1:2], data$weights)
  out[2] <- get_OAw(data[, 1:2], data$weights)
  names(out) <- c("MAEw", "OAw")
  return(out)
}

# Weighted summary for drainage classes as factor

WeightedSummary_DCfac <- function(
    data,
    lev = NULL,
    model = NULL,
    ...) {
  #Check data
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
    stop("levels of observed and predicted data do not match")
  
  require(stats)
  
  has_class_probs <- all(lev %in% colnames(data))
  
  if(has_class_probs) {
    data$pred_num <- apply(
      data[, colnames(data) %in% lev],
      1,
      function(x2) {
        out2 <- stats::weighted.mean(
          x = c(1:length(lev)),
          w = x2,
          na.rm = TRUE
          )
        return(out2)
      }
    )
  } else {
    data$pred_num <- as.numeric(pred)
  }
  get_MAEw <- function(d, w) {
    require(stats)
    if (nrow(d) == 0) {
      out <- NA
    } else {
      d$obs_num <- as.numeric(d$obs)
      ae <- abs(d$obs_num - d$pred_num)
      out <- stats::weighted.mean(x = ae, w = w, na.rm = TRUE)
    }
    return(out)
  }
  get_OAw <- function(d, w) {
    require(stats)
    if (nrow(d) == 0) {
      out <- NA
    } else {
      d_equal <- d$obs == d$pred
      out <- stats::weighted.mean(x = d_equal, w = w, na.rm = TRUE)
    }
    return(out)
  }
  out <- numeric()
  out[1] <- get_MAEw(data, data$weights)
  out[2] <- get_OAw(data, data$weights)
  names(out) <- c("MAEw", "OAw")
  return(out)
}

# END