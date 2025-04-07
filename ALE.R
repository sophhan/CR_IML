# Helper function to order feature levels
order_levels <- function(data, variable_values, variable_name) {
  feature <- droplevels(variable_values)
  x.count <- as.numeric(table(feature))
  x.prob <- x.count / sum(x.count)
  K <- nlevels(feature)
  
  dists <- lapply(setdiff(colnames(data), variable_name), function(x) {
    feature.x <- data[, x]
    dists <- expand.grid(levels(feature), levels(feature))
    colnames(dists) <- c("from.level", "to.level")
    if (inherits(feature.x, "factor")) {
      A <- table(feature, feature.x) / x.count
      dists$dist <- rowSums(abs(A[dists[, "from.level"], ] - A[dists[, "to.level"], ])) / 2
    } else {
      quants <- quantile(feature.x, probs = seq(0, 1, length.out = 100), na.rm = TRUE, names = FALSE)
      ecdfs <- data.frame(lapply(levels(feature), function(lev) {
        x.ecdf <- ecdf(feature.x[feature == lev])(quants)
      }))
      colnames(ecdfs) <- levels(feature)
      ecdf.dists.all <- abs(ecdfs[, dists$from.level] - ecdfs[, dists$to.level])
      dists$dist <- apply(ecdf.dists.all, 2, max)
    }
    dists
  })
  
  dists.cumulated.long <- Reduce(function(d1, d2) {
    d1$dist <- d1$dist + d2$dist
    d1
  }, dists)
  dists.cumulated <- xtabs(dist ~ from.level + to.level, dists.cumulated.long)
  scaled <- cmdscale(dists.cumulated, k = 1)
  order(scaled)
}

# Function to compute ALE values for all features
compute_ale <- function(model,
                        data,
                        event_type,
                        grid_points,
                        center = FALSE,
                        categorical_variables,
                        N = 500) {
  # Sample from the original data
  data <- temp[sample(1:nrow(temp), N), , drop = FALSE]
  # Convert data to data frame
  data <- as.data.frame(data)
  # Get the names of the variables
  variables <- colnames(data[, !(names(data) %in% c("status", "duration"))])
  # Convert categorical variables to factors
  data[categorical_variables] <- lapply(data[categorical_variables], factor)
  
  # Get the original predictions
  predictions_original <- predict(model, data, outcome = "test")$cif[, , event_type]
  times <- predict(model, data, outcome = "test")$time
  mean_pred <- colMeans(predictions_original)
  
  # Compute ALE profiles for each variable
  profiles <- lapply(variables, function(variable) {
    X_lower <- X_upper <- data
    variable_values <- data[, variable]
    if (variable %in% categorical_variables) {
      if (!is.factor(variable_values)) {
        is_numeric <- is.numeric(variable_values)
        is_factorized <- TRUE
        variable_values <- as.factor(variable_values)
      } else {
        is_factorized <- FALSE
      }
      levels_original <- levels(droplevels(variable_values))
      levels_n <- nlevels(droplevels(variable_values))
      
      if (inherits(variable_values, "ordered")) {
        level_order <- 1:levels_n
      } else {
        level_order <- order_levels(data, variable_values, variable)
      }
      
      # The new order of the levels
      levels_ordered <- levels_original[level_order]
      
      # The feature with the levels in the new order
      x_ordered <- order(level_order)[as.numeric(droplevels(variable_values))]
      
      # Filter rows which are not already at maximum or minimum level values
      row_ind_increase <- (1:nrow(data))[x_ordered < levels_n]
      row_ind_decrease <- (1:nrow(data))[x_ordered > 1]
      
      if (is_factorized) {
        levels_ordered <- as.character(levels_ordered)
        if (is_numeric) {
          levels_ordered <- as.numeric(levels_ordered)
        }
      }
      
      X_lower[row_ind_decrease, variable] <- levels_ordered[x_ordered[row_ind_decrease] - 1]
      X_upper[row_ind_increase, variable] <- levels_ordered[x_ordered[row_ind_increase] + 1]
      
      # Make predictions for decreased levels (excluding minimum levels)
      predictions_lower <- predict(model, X_lower[row_ind_decrease, ], outcome = "test")$cif[, , event_type]
      
      # Make predictions for increased levels (excluding maximum levels)
      predictions_upper <- predict(model, X_upper[row_ind_increase, ], outcome = "test")$cif[, , event_type]
      
      d_increase <- predictions_upper - predictions_original[row_ind_increase, ]
      d_decrease <- predictions_original[row_ind_decrease, ] - predictions_lower
      prediction_deltas <- rbind(d_increase, d_decrease)
      colnames(prediction_deltas) <- times
      
      
      deltas <- data.frame(
        interval = rep(c(x_ordered[row_ind_increase], x_ordered[row_ind_decrease] - 1), each = length(times)),
        time = rep(times, times = nrow(prediction_deltas)),
        yhat = c(t(prediction_deltas))
      )
      
      deltas <- aggregate(yhat ~ interval + time, data = deltas, FUN = mean)
      deltas1 <- deltas[deltas$interval == 1, ]
      deltas1$yhat <- 0
      deltas$interval <- deltas$interval + 1
      deltas <- rbind(deltas, deltas1)
      deltas <- deltas[order(deltas$time, deltas$interval), ]
      rownames(deltas) <- NULL
      deltas$yhat_cumsum <- ave(deltas$yhat, deltas$time, FUN = cumsum)
      
      x_count <- as.numeric(table(variable_values))
      x_prob <- x_count / sum(x_count)
      
      ale_means <- aggregate(
        yhat_cumsum ~ time,
        data = deltas,
        FUN = function(x) {
          sum(x * x_prob[level_order])
        }
      )
      colnames(ale_means)[2] <- "ale0"
      
      ale_values <- merge(deltas, ale_means, all.x = TRUE, by = "time")
      
      ale_values$ale <- ale_values$yhat_cumsum - ale_values$ale0
      ale_values$level <- levels_ordered[ale_values$interval]
      
      ale_values <- ale_values[order(ale_values$interval, ale_values$time), ]
      ale_values$ale <- ale_values$ale + mean_pred
      
      if (center) {
        min_rows <- ale_values[ale_values$level == min(ale_values$level), ]
        merged_df <- merge(ale_values,
                           min_rows[, c("time", "ale")],
                           by = "time",
                           suffixes = c("", "_min"))
        merged_df$ale <- merged_df$ale - merged_df$ale_min
        ale_values <- merged_df[, c("level", "time", "ale")]
      }
      
      return(
        data.frame(
          `_vname_` = variable,
          `_vtype_` = "categorical",
          `_x_` = ale_values$level,
          `_times_` = ale_values$time,
          `_yhat_` = ale_values$ale,
          check.names = FALSE
        )
      )
    } else {
      # Number of quantile points for determined by grid length
      quantile_vals <- as.numeric(quantile(
        variable_values,
        seq(0.01, 1, length.out = grid_points),
        type = 1
      ))
      
      # Quantile points vector
      quantile_vec <- c(min(variable_values), quantile_vals)
      quantile_vec <- unique(quantile_vec)
      quantile_df <- data.frame(id = 1:length(quantile_vec), value = quantile_vec)
      
      # Match feature instances to quantile intervals
      interval_index <- findInterval(variable_values, quantile_vec, left.open = TRUE)
      
      # Points in interval 0 should be in interval 1
      interval_index[interval_index == 0] <- 1
      
      # Prepare datasets with upper and lower interval limits replacing original feature values
      X_lower[, variable] <- quantile_vec[interval_index]
      X_upper[, variable] <- quantile_vec[interval_index + 1]
      # Get survival predictions for instances of upper and lower interval limits
      predictions_lower <- predict(model, X_lower, outcome = "test")$cif[, , event_type]
      
      predictions_upper <- predict(model, X_upper, outcome = "test")$cif[, , event_type]
      
      # First order finite differences
      prediction_deltas <- predictions_upper - predictions_lower
      # Rename columns to time points for which predictions were made
      colnames(prediction_deltas) <- times
      
      deltas <- data.frame(
        x = rep(X_lower[, variable], each = length(times)),
        interval = rep(interval_index, each = length(times)),
        time = rep(times, times = nrow(data)),
        yhat = c(t(prediction_deltas))
      )
      
      deltas <- aggregate(yhat ~ interval + time, data = deltas, FUN = mean)
      deltas$yhat_cumsum <- ave(deltas$yhat, deltas$time, FUN = cumsum)
      interval_n <- as.numeric(table(interval_index))
      n <- sum(interval_n)
      
      ale_means <- aggregate(
        yhat_cumsum ~ time,
        data = deltas,
        FUN = function(x) {
          sum(((c(0, x[1:(length(x) - 1)]) + x) / 2) * interval_n / n)
        }
      )
      colnames(ale_means)[2] <- "ale0"
      
      # Centering the ALEs to obtain final ALE values
      ale_values <- merge(deltas, ale_means, all.x = TRUE, by = "time")
      
      ale_values$ale <- ale_values$yhat_cumsum - ale_values$ale0
      ale_values$interval <- ale_values$interval + 1
      ale_values1 <- ale_values[seq(1, nrow(ale_values), length(quantile_vec) - 1), ]
      ale_values1$interval <- 1
      ale_values <- rbind(ale_values, ale_values1)
      
      ale_values <- merge(ale_values,
                          quantile_df,
                          by.x = "interval",
                          by.y = "id")
      ale_values <- ale_values[order(ale_values$interval, ale_values$time), ]
      ale_values$ale <- ale_values$ale + mean_pred

      if (center) {
        min_rows <- ale_values[ale_values$value == min(ale_values$value), ]
        merged_df <- merge(ale_values,
                           min_rows[, c("time", "ale")],
                           by = "time",
                           suffixes = c("", "_min"))
        merged_df$ale <- merged_df$ale - merged_df$ale_min
        ale_values <- merged_df[, c("value", "time", "ale")]
      }
      
      return(
        data.frame(
          `_vname_` = variable,
          `_vtype_` = "numerical",
          `_x_` = ale_values$value,
          `_times_` = ale_values$time,
          `_yhat_` = ale_values$ale,
          check.names = FALSE
        )
      )
    }
  })
  
  profiles <- do.call(rbind, profiles)
  return(profiles)
}

# Plotting function
plot_ale <- function(profiles, variable, categorical_variables) {
  plot_df <- profiles[(profiles$`_vname_` == variable),]
  plot_df$`_times_` <- as.numeric(plot_df$`_times_`)
  if (variable %in% categorical_variables) {
    # Categorical Plot
    ggplot(plot_df) +
      geom_line(aes(x = `_times_`, y = `_yhat_`, color = `_x_`, group = `_x_`)) +
      labs(title = paste0("ALE plot for ", variable),
           x = "time",
           y = "prediction") +
      theme_minimal() +
      guides(color=guide_legend(title="feature value"))
  } else {
    # Numerical Plot 
    plot_df$`_x_` <- as.numeric(plot_df$`_x_`)
    ggplot(plot_df) +
      geom_line(aes(x = `_x_`, y = `_yhat_`, color = `_times_`, group = `_times_`)) +
      scale_color_gradient() +
      labs(title = paste0("ALE plot for ", variable),
           x = variable,
           y = "prediction") +
      theme_minimal() +
      guides(color = guide_colorbar(title = "time"))
  }
}


# Compute ALE for the credit risk model
categorical_variables <- c(
  "CAR_TYPE",
  "BRANCH_services",
  "CUSTOMER_BEFORE",
  "BRANCH_health",
  "firm_LTD_PARTN",
  "BRANCH_other",
  "BRANCH_Sale",
  "PHONE",
  "REGION3",
  "BRANCH_Production",
  "REGION1",
  "REGION2",
  "BRANCH_construction"
)
# default
# ale_default <- compute_ale(
#   model = v.obj5,
#   data = temp,
#   event_type = 1,
#   grid_points = 20,
#   categorical_variables,
#   center = FALSE,
#   N = 1000
# )
ale_default_c <- compute_ale(
  model = v.obj5,
  data = temp,
  event_type = 1,
  grid_points = 20,
  categorical_variables,
  center = TRUE,
  N = 1000
)
# prepayment
# ale_prepayment <- compute_ale(
#   model = v.obj5,
#   temp,
#   event_type = 2,
#   grid_points = 20,
#   categorical_variables,
#   center = FALSE,
#   N = 1000
# )
ale_prepayment_c <- compute_ale(
  model = v.obj5,
  temp,
  event_type = 2,
  grid_points = 20,
  categorical_variables,
  center = TRUE,
  N = 1000
)

# plot ale values default
plot_default_c <- cowplot::plot_grid(
  plot_ale(ale_default_c, "NET_ASSETS", categorical_variables),
  plot_ale(ale_default_c, "DOWNPAYMENT_TO_CAR_PRICE", categorical_variables),
  plot_ale(ale_default_c, "BRANCH_services", categorical_variables),
  plot_ale(ale_default_c, "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", categorical_variables),
  plot_ale(ale_default_c, "COMPANY_AGE", categorical_variables),
  plot_ale(ale_default_c, "CAR_PRICE", categorical_variables),
  plot_ale(ale_default_c, "CUSTOMER_BEFORE", categorical_variables),
  plot_ale(ale_default_c, "ANNUAL_TURNOVER_LAST_YEAR", categorical_variables),
  plot_ale(ale_default_c, "BRANCH_other", categorical_variables),
  nrow = 3, ncol = 3)
plot_default_c
ggsave(
  "/Users/sophielangbein/Desktop/PhD/CR_IML/plots/ale_default_c.pdf",
  plot = plot_default_c,
  width = 15,
  height = 10
)
# plot centered ale values default
plot_prepayment_c <- cowplot::plot_grid(
  plot_ale(ale_prepayment_c, "NET_ASSETS", categorical_variables),
  plot_ale(ale_prepayment_c, "car_age", categorical_variables),
  plot_ale(ale_prepayment_c, "COMPANY_AGE", categorical_variables),
  plot_ale(ale_prepayment_c, "PHONE", categorical_variables),
  plot_ale(ale_prepayment_c, "ANNUAL_TURNOVER_LAST_YEAR", categorical_variables),
  plot_ale(ale_prepayment_c, "ANNUAL_INCOME_LAST_YEAR", categorical_variables),
  plot_ale(ale_prepayment_c, "NUMBER_OF_EMPLOYEES", categorical_variables),
  plot_ale(ale_prepayment_c, "length_coop", categorical_variables),
  plot_ale(ale_prepayment_c, "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", categorical_variables),
  nrow = 3, ncol = 3)
plot_prepayment_c
ggsave(
  "/Users/sophielangbein/Desktop/PhD/CR_IML/plots/ale_prepayment_c.pdf",
  plot = plot_prepayment_c,
  width = 15,
  height = 10
)
plot_ale(ale_prepayment_c, "PHONE", categorical_variables)
