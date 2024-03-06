# Informations ------------------------------------------------------------
# Title: reox_model.R
# Author: CHAMOUX Maxime, SOLSONNA Robert & BOUDRY Félix
# License: GPLv3
# Description: Vastus lateralis reoxygenation model.

set.seed(123) # Using fixed seed for reproducibility

# Loading required packages
library(nlstools)
library(broom)
library(ggplot2)
library(tidyverse)
library(nls2)
library(readxl)
library(rlist)
library(ggpubr)
library(parallel)

# Select parameters
# Choose those parameters carefully regarding your hardware configuration,
# nls is using HUGE amounts of memory
n_trials <- 300000 # Number of models to try
n_cores <- 5 # number of cores to use
data_path <- "Data/BFR_reox"

# Model's formula
nls_formula <-
  tsi ~ tsi_base +
  ifelse(
    test = time <= extracted_value,
    yes = 0,
    no = A1 * (1 - exp(-(
      time - extracted_value
    ) / tau1))
  ) +
  ifelse(test = time <= TD2,
         yes = 0,
         no = A2 * (1 - exp(-(time - TD2) / tau2)))

# Modeling
nls_params <- data.frame(
  tau1 = c(0, 60),
  tau2 = c(0, 200),
  A1 = c(0, 80),
  A2 = c(-50, 50),
  TD2 = c(80, 250)
)

# Functions
nls_fit <- function(df_name) {
  target_df <- reox_data[[df_name]]

  # Select starting value
  extracted_value <<- model_start$td1[model_start$name == df_name]

  tsi_base <<-
    target_df %>%
    slice(which(time == floor(extracted_value))) %>%
    select(tsi) %>%
    as.numeric()

  # Fit model
  result_model <- nls2(
    formula = nls_formula,
    data = target_df,
    start = nls_params,
    algorithm = "random-search",
    control = nls.control(maxiter = n_trials)
  )
  return(result_model)
}

# Import and format data
reox_data_path <-
  list.files(
    path = data_path,
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  )
reox_data_name <-
  list.files(
    path = data_path,
    pattern = "\\.xlsx",
    all.files = TRUE,
    full.names = TRUE
  ) %>%
  tools::file_path_sans_ext() %>%
  basename()
reox_data <- lapply(reox_data_path, read_xlsx) %>%
  `names<-`(x = .,
            value = reox_data_name) %>%
  lapply(setNames, c("time", "tsi"))
model_start <- read_xlsx("Data/td1_BFR.xlsx") %>%
  as.data.frame()

reox_data <- lapply(names(reox_data), function(df_name) {
  # Find corresponding start value in 'model_start'
  target_df <- reox_data[[df_name]]
  start_value <- model_start$td1[model_start$name == df_name] - 1
  filtered_data <- filter(target_df, time >= floor(start_value))
  return(filtered_data)
}) %>%
  `names<-`(value = reox_data_name)

nls_results <- lapply(names(reox_data), nls_fit)

# Results
TSS <-
  as.data.frame(list.rbind(lapply(reox_data, function(target_df) {
    sum((target_df$tsi - mean(target_df$tsi)) ^ 2)
  }))) %>%
  `colnames<-`(value = "TSS")

RSS <-
  as.data.frame(list.rbind(lapply(nls_results, function(target_df) {
    residuals <- residuals(target_df)
    RSS <- sum(residuals ^ 2)
  }))) %>%
  `colnames<-`(value = "RSS")

rsquared <- cbind(RSS, TSS) %>%
  mutate(rsquared = 1 - (RSS / TSS))
predicted_params <-
  as.data.frame(list.rbind(lapply(nls_results, function(target_df) {
    coef(target_df)
  })))

nls_metrics <- as.data.frame(cbind(predicted_params, rsquared))

# Plot
predicted_results <- lapply(nls_results, function(target_df) {
  augment(target_df)
}) %>%
  `names<-`(value = reox_data_name)

plot_list <- mapply(
  function(model_results, raw_data, model_metrics) {
    plot <- ggplot(data = raw_data) +
      geom_point(aes(time, tsi)) +
      geom_line(
        aes(time, .fitted),
        data = model_results,
        color = "red",
        linewidth = 1.5
      ) +
      annotate(
        "text",
        x = 180,
        y = -Inf,
        label = paste("R-squared:", round(model_metrics, 3)),
        color = "black",
        size = 3,
        vjust = -0.5
      )
    return(plot)
  },
  model_results = predicted_results,
  raw_data = reox_data,
  model_metrics = nls_metrics$rsquared,
  SIMPLIFY = FALSE
)

result_plot <-
  ggarrange(
    plotlist = plot_list,
    labels = reox_data_name,
    font.label = list(size = 8, color = "black"),
    vjust = 0.2
  ) +
  theme(plot.margin = margin(t = 0.5, r = 0.1, b = 0.5, l = 0.1, "cm"))
