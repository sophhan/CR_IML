# Function to compute PDP values for all features
compute_pdp <- function(model,
                        data,
                        event_type,
                        grid_points,
                        center = FALSE,
                        categorical_variables,
                        N = 500) {
  # Sample data
  data = as.data.frame(data)
  data = data[sample(nrow(data), N), ]
  # Initialize values
  variables <- colnames(data[, !(names(data) %in% c("status", "duration"))])
  variable_splits_type = "quantiles"
  times <- predict(model, data, outcome = "test")$time
  
  # Function to calculate variable split values 
  calculate_variable_split <- function(data,
                                       variables = colnames(data),
                                       categorical_variables = NULL,
                                       grid_points = 101,
                                       variable_splits_type = "quantiles",
                                       new_observation = NA) {
    variable_splits <- lapply(variables, function(var) {
      selected_column <- na.omit(data[, var])
      
      if (!(var %in% categorical_variables)) {
        probs <- seq(0, 1, length.out = grid_points)
        if (variable_splits_type == "quantiles") {
          selected_splits <- unique(quantile(selected_column, probs = probs))
        } else {
          selected_splits <- seq(
            min(selected_column, na.rm = TRUE),
            max(selected_column, na.rm = TRUE),
            length.out = grid_points
          )
        }
        if (!any(is.na(new_observation))) {
          selected_splits <- sort(unique(c(
            selected_splits, na.omit(new_observation[, var])
          )))
        }
      } else {
        if (any(is.na(new_observation))) {
          selected_splits <- sort(unique(selected_column))
        } else {
          selected_splits <- sort(unique(rbind(data[, var, drop = FALSE], new_observation[, var, drop = FALSE])[, 1]))
        }
      }
      selected_splits
    })
    names(variable_splits) <- variables
    variable_splits
  }
  
  # Calculate variable splits
  variable_splits <- calculate_variable_split(
    data,
    variables = variables,
    categorical_variables = categorical_variables,
    grid_points = grid_points,
    variable_splits_type = variable_splits_type,
    new_observation = NA
  )
  
  # Name variable splits
  variables <- names(variable_splits)
  
  # Calculate original predictions and times
  predictions_original <- predict(model, data, outcome = "test")$cif[, , event_type]
  times <- predict(model, data, outcome = "test")$time
  
  # Compute PDP values for each variable
  profiles <- lapply(variables, function(variable) {
    split_points <- variable_splits[[variable]]
    
    new_data <- data[rep(1:nrow(data), each = length(split_points)), , drop = FALSE]
    new_data[, variable] <- rep(split_points, nrow(data))
    new_data_sorted <- new_data[order(new_data[,variable]), ]
    
    yhat = predict(model, new_data_sorted, outcome = "test")$cif[, , event_type]
    group_index <- rep(1:length(split_points), each = N)
    yhat_avg <- apply(yhat, 2, function(col) tapply(col, group_index, mean))
    if (center) {
      yhat_avg <- sweep(yhat_avg, 2, yhat_avg[1, ], FUN = "-")
    }
    result <- data.frame(
      `_times_` = rep(times, times = length(split_points)),
      `_x_` = rep(split_points, each = length(times)),
      `_yhat_` = c(t(yhat_avg)),
      `_vname_` = variable,
      check.names = FALSE
    )
    result
  })
  
  # Combine profiles
  profile <- do.call(rbind, profiles)
  
  return(profile)
}

# List of categorical variables
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

# Plotting function
plot_pdp <- function(profiles, variable, categorical_variables) {
  plot_df <- profiles[(profiles$`_vname_` == variable),]
  plot_df$`_times_` <- as.numeric(plot_df$`_times_`)
  if (variable %in% categorical_variables) {
    plot_df[, "_x_"] <- as.factor(plot_df[, "_x_"])
    # Categorical Plot
    ggplot(plot_df) +
      geom_line(aes(x = `_times_`, y = `_yhat_`, color = `_x_`, group = `_x_`)) +
      labs(title = paste0("PDP for ", variable),
           x = "time",
           y = "prediction") +
      theme_minimal() +
      scale_color_discrete() +
      guides(color=guide_legend(title="feature value"))
  } else {
    # Numerical Plot 
    plot_df$`_x_` <- as.numeric(plot_df$`_x_`)
    ggplot(plot_df) +
      geom_line(aes(x = `_x_`, y = `_yhat_`, color = `_times_`, group = `_times_`)) +
      scale_color_gradient() +
      labs(title = paste0("PDP for ", variable),
           x = variable,
           y = "prediction") +
      theme_minimal() +
      guides(color = guide_colorbar(title = "time"))
  }
}


# default
# pdp_default <- compute_pdp(
#   model = v.obj5,
#   data = as.data.frame(temp),
#   event_type = 1,
#   grid_points = 20,
#   categorical_variables,
#   center = FALSE,
#   N = 1000
# )
pdp_default_c <- compute_pdp(
  model = v.obj5,
  data = as.data.frame(temp),
  event_type = 1,
  grid_points = 20,
  categorical_variables,
  center = TRUE,
  N = 1000
)

# prepayment
# pdp_prepayment <- compute_pdp(
#   model = v.obj5,
#   data = as.data.frame(temp),
#   event_type = 2,
#   grid_points = 20,
#   categorical_variables,
#   center = FALSE,
#   N = 500
# )
pdp_prepayment_c <- compute_pdp(
  model = v.obj5,
  data = as.data.frame(temp),
  event_type = 2,
  grid_points = 20,
  categorical_variables,
  center = TRUE,
  N = 1000
)


# plot pdp values default
# plot_default <- cowplot::plot_grid(
#   plot_pdp(pdp_default, "CAR_TYPE", categorical_variables),
#   plot_pdp(pdp_default, "NET_ASSETS", categorical_variables),
#   plot_pdp(pdp_default, "DOWNPAYMENT_TO_CAR_PRICE", categorical_variables),
#   nrow = 3)
# plot_default

# plot pdp values prepayment 
# plot_prepayment <- cowplot::plot_grid(
#   plot_pdp(pdp_prepayment, "CAR_TYPE", categorical_variables),
#   plot_pdp(pdp_prepayment, "NET_ASSETS", categorical_variables),
#   plot_pdp(pdp_prepayment, "DOWNPAYMENT_TO_CAR_PRICE", categorical_variables),
#   nrow = 3)
# plot_prepayment

# plot pdp values default
# plot_default <- cowplot::plot_grid(
#   plot_pdp(pdp_default, "NET_ASSETS", categorical_variables),
#   plot_pdp(pdp_default, "DOWNPAYMENT_TO_CAR_PRICE", categorical_variables),
#   plot_pdp(pdp_default, "BRANCH_services", categorical_variables),
#   plot_pdp(pdp_default, "ANNUAL_TURNOVER_LAST_YEAR", categorical_variables),
#   plot_pdp(pdp_default, "NUMBER_OF_EMPLOYEES", categorical_variables),
#   plot_pdp(pdp_default, "COMPANY_AGE", categorical_variables),
#   plot_pdp(pdp_default, "ANNUAL_INCOME_LAST_YEAR", categorical_variables),
#   plot_pdp(pdp_default, "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", categorical_variables),
#   plot_pdp(pdp_default, "BRANCH_other", categorical_variables),
#   nrow = 3, ncol = 3)
# plot_default
plot_default_c <- cowplot::plot_grid(
  plot_pdp(pdp_default_c, "NET_ASSETS", categorical_variables),
  plot_pdp(pdp_default_c, "DOWNPAYMENT_TO_CAR_PRICE", categorical_variables),
  plot_pdp(pdp_default_c, "BRANCH_services", categorical_variables),
  plot_pdp(pdp_default_c, "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", categorical_variables),
  plot_pdp(pdp_default_c, "COMPANY_AGE", categorical_variables),
  plot_pdp(pdp_default_c, "CAR_PRICE", categorical_variables),
  plot_pdp(pdp_default_c, "CUSTOMER_BEFORE", categorical_variables),
  plot_pdp(pdp_default_c, "ANNUAL_TURNOVER_LAST_YEAR", categorical_variables),
  plot_pdp(pdp_default_c, "BRANCH_other", categorical_variables),
  nrow = 3, ncol = 3)
plot_default_c
ggsave(
  "/Users/sophielangbein/Desktop/PhD/CR_IML/plots/pdp_default_c.pdf",
  plot = plot_default_c,
  width = 15,
  height = 10
)
# plot centered pdp values default
plot_prepayment_c <- cowplot::plot_grid(
  plot_pdp(pdp_prepayment_c, "NET_ASSETS", categorical_variables),
  plot_pdp(pdp_prepayment_c, "car_age", categorical_variables),
  plot_pdp(pdp_prepayment_c, "COMPANY_AGE", categorical_variables),
  plot_pdp(pdp_prepayment_c, "PHONE", categorical_variables),
  plot_pdp(pdp_prepayment_c, "ANNUAL_TURNOVER_LAST_YEAR", categorical_variables),
  plot_pdp(pdp_prepayment_c, "ANNUAL_INCOME_LAST_YEAR", categorical_variables),
  plot_pdp(pdp_prepayment_c, "NUMBER_OF_EMPLOYEES", categorical_variables),
  plot_pdp(pdp_prepayment_c, "length_coop", categorical_variables),
  plot_pdp(pdp_prepayment_c, "NUMBER_OF_ALL_BAD_PREVIOUS_CONTR", categorical_variables),
  nrow = 3, ncol = 3)
plot_prepayment_c
ggsave(
  "/Users/sophielangbein/Desktop/PhD/CR_IML/plots/pdp_prepayment_c.pdf",
  plot = plot_prepayment_c,
  width = 15,
  height = 10
)

