# Define a function to compute the cause-specific Brier Score
compute_brier_score <- function(model, data, event_type) {
  # Step 1: Compute the Kaplan-Meier estimate for censoring (IPCW)
  km_censor <- survfit(Surv(duration, status == 0) ~ 1, data = data)
  
  # Function to get IPCW weights
  get_ipcw <- function(t) {
    censor_prob <- summary(km_censor, times = t, extend = TRUE)$surv
    ifelse(censor_prob == 0, NA, 1 / censor_prob)  # Avoid division by zero
  }
  
  # Step 2: Compute predicted cumulative incidence for event_type
  times <- predict(model, data, outcome = "test")$time
  pred_cif <- predict(model, data, outcome = "test")$cif[, , event_type]
  
  # Step 3: Initialize Brier score storage
  brier_scores <- numeric(length(times))
  
  for (j in seq_along(times)) {
    t <- times[j]
    
    # Step 4: Compute the indicator function for observed events
    event_observed <- (data$duration <= t) &
      (data$status == event_type)
    
    # Step 5: Get IPCW weights for censoring
    weights <- get_ipcw(t)
    
    # Step 6: Compute the Brier Score at time t
    brier_scores[j] <- mean(weights * (pred_cif[, j] - event_observed)^2, na.rm = TRUE)
  }
  
  return(data.frame(time = times, brier = brier_scores))
}

# Define a function to compute pfi using the cause-specific Brier Score
compute_pfi <- function(model, data, event_type) {
  brier_score_full <- compute_brier_score(model, data, event_type)
  col_names <- colnames(data[, setdiff(colnames(data), c("status", "duration"))])
  data_df <- as.data.frame(data)
  
  # Function to compute PFI for a single feature
  pfi <- function(list_element) {
    data_perm <- data_df
    data_perm[, list_element] <- sample(data_perm[, list_element])
    model_perm <- rfsrc(
      Surv(duration, status) ~ .,
      data = data_perm,
      ntree = 1000,
      nodesize = 15,
      mtry = 6,
      nsplit = 10,
      samptype = "swr",
      splitrule = "logrankCR",
      importance = TRUE,
      statistics = TRUE
    )
    brier_score_perm <- compute_brier_score(model_perm, data_perm, event_type)
    return(brier_score_perm[, "brier"] - brier_score_full[, "brier"])
  }
  
  # Compute PFI for all features
  pfi_list <- lapply(col_names, pfi)
  times <- brier_score_full$time
  
  # Convert to long format dataframe
  pfi_df <- do.call(rbind, lapply(seq_along(pfi_list), function(i) {
    data.frame(pfi_value = pfi_list[[i]],
               feature = col_names[i],
               time = times)
  }))
  
  return(pfi_df)
}


### Compute PFI for rfsrc credit risk model
# default
pfi_cr_event1 <- compute_pfi(model = v.obj5, data = temp, event_type = 1)
# prepayment
pfi_cr_event2 <- compute_pfi(model = v.obj5, data = temp, event_type = 2)


### Plot PFI
# default
# Compute the average pfi_value over time
plot_pfi <- function(data, n_features = NULL, title = "Default Prediction") {
  if (!is.null(n_features)) {
    df_avg <- data %>%
      group_by(feature) %>%
      summarise(avg_pfi_value = mean(pfi_value, na.rm = TRUE)) %>%
      arrange(desc(avg_pfi_value))
    # Select the top n features
    topn <- df_avg$feature[1:n_features]
    # Filter the dataframe to include only the top n features
    data <- data %>%
      filter(feature %in% topn)
  }
  # Plot the pfi_values over time
  ggplot(data,
         aes(
           x = time,
           y = pfi_value,
           color = feature,
           linetype = feature
         )) +
    geom_line() +
    labs(title = paste("Permutation Feature Importance for", title) , 
         x = "Time", 
         y = "Brier Score Difference") +
    theme_minimal()
}

# Plot results
pfi_cr_default_10 <- plot_pfi(pfi_cr_event1, n_features = 10)
pfi_cr_default_10
ggsave(
  "/Users/sophielangbein/Desktop/PhD/CR_IML/plots/pfi_default_10.pdf",
  plot = pfi_cr_default_10,
  width = 10,
  height = 4
)

pfi_cr_prepayment_10 <- plot_pfi(pfi_cr_event2, n_features = 10, title = "Prepayment Prediction")
pfi_cr_prepayment_10
ggsave(
  "/Users/sophielangbein/Desktop/PhD/CR_IML/plots/pfi_prepayment_10.pdf",
  plot = pfi_cr_prepayment_10,
  width = 10,
  height = 4
)

