# Residential Construction Time Series Analysis
# This script benchmarks multiple Holt-Winters specifications using a strict
# out-of-sample backtest and exports artifacts for auditability.

options(stringsAsFactors = FALSE, scipen = 999)

# Package validation is explicit to keep runtime behavior deterministic across machines.
required_packages <- c(
  "readr",
  "dplyr",
  "lubridate",
  "forecast",
  "Metrics",
  "ggplot2"
)

missing_packages <- setdiff(required_packages, rownames(utils::installed.packages()))
if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Missing packages detected: %s. Install them before running the analysis.",
      paste(missing_packages, collapse = ", ")
    )
  )
}

get_script_directory <- function() {
  script_flag <- "--file="
  matching_args <- grep(script_flag, commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(matching_args) == 0L) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  script_path <- sub(script_flag, "", matching_args[1], fixed = TRUE)
  dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
}

parse_monthly_date <- function(month_string) {
  current_locale <- Sys.getlocale("LC_TIME")
  on.exit(Sys.setlocale("LC_TIME", current_locale), add = TRUE)
  Sys.setlocale("LC_TIME", "C")
  as.Date(paste0("01-", month_string), format = "%d-%b-%y")
}

compute_error_metrics <- function(actual, predicted) {
  actual <- as.numeric(actual)
  predicted <- as.numeric(predicted)

  if (length(actual) != length(predicted)) {
    stop("Metric computation failed: actual and predicted vectors have different lengths.")
  }

  valid_rows <- is.finite(actual) & is.finite(predicted)
  actual <- actual[valid_rows]
  predicted <- predicted[valid_rows]

  if (length(actual) == 0L) {
    stop("Metric computation failed: no valid observations after NA filtering.")
  }

  epsilon <- .Machine$double.eps
  mape <- mean(abs((actual - predicted) / pmax(abs(actual), epsilon))) * 100
  smape <- mean((2 * abs(actual - predicted)) / pmax(abs(actual) + abs(predicted), epsilon)) * 100

  data.frame(
    MAE = Metrics::mae(actual, predicted),
    MSE = Metrics::mse(actual, predicted),
    RMSE = sqrt(Metrics::mse(actual, predicted)),
    MAPE = mape,
    sMAPE = smape,
    check.names = FALSE
  )
}

fit_holt_winters <- function(series, seasonal_type, alpha = NA_real_, beta = NA_real_, gamma = NA_real_) {
  if (is.na(alpha) && is.na(beta) && is.na(gamma)) {
    return(stats::HoltWinters(series, seasonal = seasonal_type))
  }
  stats::HoltWinters(
    series,
    seasonal = seasonal_type,
    alpha = alpha,
    beta = beta,
    gamma = gamma
  )
}

build_monthly_ts <- function(data_frame) {
  stats::ts(
    data_frame$nstarted,
    start = c(lubridate::year(data_frame$date[1]), lubridate::month(data_frame$date[1])),
    frequency = 12
  )
}

# ---- Data loading and quality controls ----

project_dir <- get_script_directory()
data_path <- file.path(project_dir, "construction_data.csv")

if (!file.exists(data_path)) {
  stop(sprintf("Input file not found: %s", data_path))
}

construction_df <- readr::read_csv(
  data_path,
  col_types = readr::cols(
    date = readr::col_character(),
    nstarted = readr::col_double()
  )
)

construction_df <- dplyr::mutate(construction_df, date = parse_monthly_date(date))
construction_df <- dplyr::arrange(construction_df, date)

if (anyNA(construction_df$date)) {
  stop("Date parsing failed for at least one record.")
}
if (anyNA(construction_df$nstarted)) {
  stop("Missing target values detected in nstarted.")
}
if (anyDuplicated(construction_df$date) > 0L) {
  stop("Duplicate monthly timestamps detected. Each month must be unique.")
}
expected_months <- seq.Date(min(construction_df$date), max(construction_df$date), by = "month")
missing_months <- setdiff(expected_months, construction_df$date)
if (length(missing_months) > 0L) {
  stop(sprintf("Time index contains %d missing month(s).", length(missing_months)))
}
if (any(construction_df$nstarted <= 0)) {
  stop("Non-positive values detected. Multiplicative Holt-Winters requires strictly positive series.")
}

# We enforce a fixed holdout horizon to keep model comparison stable and auditable.
forecast_horizon <- 12L
history_window_months <- 24L

if (nrow(construction_df) <= forecast_horizon + history_window_months) {
  stop("Insufficient data to run a 12-month backtest with a 24-month history window.")
}

split_index <- nrow(construction_df) - forecast_horizon
train_df <- construction_df[seq_len(split_index), , drop = FALSE]
test_df <- construction_df[(split_index + 1L):nrow(construction_df), , drop = FALSE]

# Backtest leakage guards: holdout must be strictly out-of-time and non-overlapping.
if (nrow(test_df) != forecast_horizon) {
  stop("Backtest split error: test window length does not match forecast horizon.")
}
if (max(train_df$date) >= min(test_df$date)) {
  stop("Backtest split error: train and test periods are not strictly time-ordered.")
}
if (length(intersect(train_df$date, test_df$date)) > 0L) {
  stop("Backtest split error: train and test date ranges overlap.")
}

train_ts <- build_monthly_ts(train_df)
test_ts <- build_monthly_ts(test_df)
full_ts <- build_monthly_ts(construction_df)

# ---- Model specification and estimation ----

# A tabular model registry avoids repetitive code and simplifies governance reviews.
model_specs <- data.frame(
  model_key = c(
    "HW_Add_Auto",
    "HW_Add_Fixed",
    "HW_Mult_Auto",
    "HW_Mult_Fixed",
    "Naive",
    "Seasonal_Naive"
  ),
  display_name = c(
    "HW Additive (Auto)",
    "HW Additive (Fixed)",
    "HW Multiplicative (Auto)",
    "HW Multiplicative (Fixed)",
    "Naive",
    "Seasonal Naive"
  ),
  model_family = c(
    "holt_winters",
    "holt_winters",
    "holt_winters",
    "holt_winters",
    "naive",
    "snaive"
  ),
  seasonal_type = c("additive", "additive", "multiplicative", "multiplicative", NA, NA),
  alpha = c(NA_real_, 0.3, NA_real_, 0.4, NA_real_, NA_real_),
  beta = c(NA_real_, 0.1, NA_real_, 0.2, NA_real_, NA_real_),
  gamma = c(NA_real_, 0.4, NA_real_, 0.3, NA_real_, NA_real_),
  stringsAsFactors = FALSE
)

model_results <- lapply(seq_len(nrow(model_specs)), function(i) {
  spec <- model_specs[i, , drop = FALSE]
  if (spec$model_family == "holt_winters") {
    fitted_model <- fit_holt_winters(
      series = train_ts,
      seasonal_type = spec$seasonal_type,
      alpha = spec$alpha,
      beta = spec$beta,
      gamma = spec$gamma
    )
    model_forecast <- forecast::forecast(fitted_model, h = forecast_horizon)
    hw_fitted <- fitted_model$fitted[, "xhat"]
    fitted_values <- as.numeric(hw_fitted)
    aligned_train_actual <- as.numeric(
      stats::window(train_ts, start = stats::start(hw_fitted), end = stats::end(hw_fitted))
    )

    estimated_alpha <- ifelse(is.na(spec$alpha), fitted_model$alpha, spec$alpha)
    estimated_beta <- ifelse(is.na(spec$beta), fitted_model$beta, spec$beta)
    estimated_gamma <- ifelse(is.na(spec$gamma), fitted_model$gamma, spec$gamma)
  } else if (spec$model_family == "naive") {
    fitted_model <- forecast::naive(train_ts, h = forecast_horizon)
    model_forecast <- fitted_model
    fitted_values <- as.numeric(fitted_model$fitted)
    aligned_train_actual <- as.numeric(train_ts)

    estimated_alpha <- NA_real_
    estimated_beta <- NA_real_
    estimated_gamma <- NA_real_
  } else if (spec$model_family == "snaive") {
    fitted_model <- forecast::snaive(train_ts, h = forecast_horizon)
    model_forecast <- fitted_model
    fitted_values <- as.numeric(fitted_model$fitted)
    aligned_train_actual <- as.numeric(train_ts)

    estimated_alpha <- NA_real_
    estimated_beta <- NA_real_
    estimated_gamma <- NA_real_
  } else {
    stop(sprintf("Unsupported model family: %s", spec$model_family))
  }

  train_metrics <- compute_error_metrics(aligned_train_actual, fitted_values)
  test_metrics <- compute_error_metrics(test_ts, model_forecast$mean)

  residual_vector <- aligned_train_actual - fitted_values
  residual_vector <- residual_vector[is.finite(residual_vector)]
  ljung_box_p_value <- NA_real_
  if (length(residual_vector) > 12L) {
    ljung_box_p_value <- stats::Box.test(residual_vector, lag = 12, type = "Ljung-Box")$p.value
  }

  list(
    spec = spec,
    model = fitted_model,
    forecast = model_forecast,
    train_metrics = train_metrics,
    test_metrics = test_metrics,
    alpha = estimated_alpha,
    beta = estimated_beta,
    gamma = estimated_gamma,
    diagnostics = data.frame(
      residual_mean = mean(residual_vector),
      residual_sd = stats::sd(residual_vector),
      ljung_box_p_value = ljung_box_p_value,
      check.names = FALSE
    )
  )
})

names(model_results) <- model_specs$model_key

build_metrics_row <- function(model_key, result_object, dataset_label, metric_frame) {
  spec <- result_object$spec
  parameter_source <- if (
    spec$model_family == "holt_winters" && is.na(spec$alpha)
  ) {
    "auto"
  } else if (spec$model_family == "holt_winters") {
    "fixed"
  } else {
    "baseline"
  }

  data.frame(
    model = model_key,
    model_name = spec$display_name,
    model_family = spec$model_family,
    seasonal_type = spec$seasonal_type,
    parameter_source = parameter_source,
    dataset = dataset_label,
    alpha = result_object$alpha,
    beta = result_object$beta,
    gamma = result_object$gamma,
    MAE = metric_frame$MAE,
    MSE = metric_frame$MSE,
    RMSE = metric_frame$RMSE,
    MAPE = metric_frame$MAPE,
    sMAPE = metric_frame$sMAPE,
    check.names = FALSE
  )
}

metrics_table <- do.call(
  rbind,
  lapply(names(model_results), function(model_key) {
    result_object <- model_results[[model_key]]
    rbind(
      build_metrics_row(model_key, result_object, "train", result_object$train_metrics),
      build_metrics_row(model_key, result_object, "test", result_object$test_metrics)
    )
  })
)

metrics_table <- dplyr::arrange(metrics_table, dataset, MAE)

cat("\nModel Performance Summary:\n")
print(metrics_table, row.names = FALSE)

diagnostics_table <- do.call(
  rbind,
  lapply(names(model_results), function(model_key) {
    result_object <- model_results[[model_key]]
    diagnostics <- result_object$diagnostics
    data.frame(
      model = model_key,
      model_name = result_object$spec$display_name,
      residual_mean = diagnostics$residual_mean,
      residual_sd = diagnostics$residual_sd,
      ljung_box_p_value = diagnostics$ljung_box_p_value,
      check.names = FALSE
    )
  })
)

cat("\nResidual Diagnostics (train window):\n")
print(diagnostics_table, row.names = FALSE)

test_metrics_table <- metrics_table[metrics_table$dataset == "test", , drop = FALSE]
champion_model_key <- test_metrics_table$model[which.min(test_metrics_table$MAE)]
cat(sprintf("\nChampion model based on test MAE: %s\n", champion_model_key))

champion_spec <- model_specs[model_specs$model_key == champion_model_key, , drop = FALSE]
if (champion_spec$model_family == "holt_winters") {
  champion_model_full <- fit_holt_winters(
    series = full_ts,
    seasonal_type = champion_spec$seasonal_type,
    alpha = champion_spec$alpha,
    beta = champion_spec$beta,
    gamma = champion_spec$gamma
  )
  champion_forecast <- forecast::forecast(champion_model_full, h = forecast_horizon)
} else if (champion_spec$model_family == "naive") {
  champion_model_full <- forecast::naive(full_ts, h = forecast_horizon)
  champion_forecast <- champion_model_full
} else if (champion_spec$model_family == "snaive") {
  champion_model_full <- forecast::snaive(full_ts, h = forecast_horizon)
  champion_forecast <- champion_model_full
} else {
  stop(sprintf("Unsupported champion model family: %s", champion_spec$model_family))
}

future_dates <- seq.Date(
  from = seq.Date(max(construction_df$date), by = "month", length.out = 2L)[2],
  by = "month",
  length.out = forecast_horizon
)

forward_forecast_table <- data.frame(
  date = future_dates,
  point_forecast = as.numeric(champion_forecast$mean),
  lower_80 = as.numeric(champion_forecast$lower[, "80%"]),
  upper_80 = as.numeric(champion_forecast$upper[, "80%"]),
  lower_95 = as.numeric(champion_forecast$lower[, "95%"]),
  upper_95 = as.numeric(champion_forecast$upper[, "95%"])
)

cat("\nChampion Model Forward Forecast (12 Months):\n")
print(forward_forecast_table, row.names = FALSE)

# ---- Visualization ----

historical_plot <- ggplot2::ggplot(construction_df, ggplot2::aes(x = date, y = nstarted)) +
  ggplot2::geom_line(color = "#1B4965", linewidth = 1) +
  ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
  ggplot2::labs(
    title = "Monthly Residential Construction Starts",
    subtitle = "Raw historical observations",
    x = "Date",
    y = "Construction starts (thousand units)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

history_df <- tail(train_df, history_window_months)
actual_frame <- dplyr::bind_rows(
  dplyr::transmute(history_df, date = date, series = "Actual", value = nstarted),
  dplyr::transmute(test_df, date = date, series = "Actual", value = nstarted)
)

forecast_frame <- dplyr::bind_rows(
  lapply(names(model_results), function(model_key) {
    data.frame(
      date = test_df$date,
      series = model_results[[model_key]]$spec$display_name,
      value = as.numeric(model_results[[model_key]]$forecast$mean),
      stringsAsFactors = FALSE
    )
  })
)

# Add one anchor point per model at the train endpoint to make forecast lines
# visually start from the last known historical observation.
forecast_anchor_frame <- data.frame(
  date = rep(max(train_df$date), times = nrow(model_specs)),
  series = model_specs$display_name,
  value = rep(tail(train_df$nstarted, 1), times = nrow(model_specs)),
  stringsAsFactors = FALSE
)

backtest_plot_data <- dplyr::bind_rows(actual_frame, forecast_anchor_frame, forecast_frame)

backtest_plot <- ggplot2::ggplot(
  backtest_plot_data,
  ggplot2::aes(x = date, y = value, color = series, linetype = series)
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_vline(
    xintercept = max(train_df$date),
    color = "#9B2226",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(
    title = "Backtest: Actuals vs Model Forecasts",
    subtitle = "Dashed line marks forecast origin (last historical observation)",
    x = "Date",
    y = "Construction starts (thousand units)",
    color = "Series",
    linetype = "Series"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

print(historical_plot)
print(backtest_plot)

# ---- Artifact export ----

output_dir <- file.path(project_dir, "analysis_outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

readr::write_csv(metrics_table, file.path(output_dir, "model_metrics.csv"))
readr::write_csv(metrics_table, file.path(output_dir, "holt_winters_model_metrics.csv"))
readr::write_csv(diagnostics_table, file.path(output_dir, "model_residual_diagnostics.csv"))
readr::write_csv(forward_forecast_table, file.path(output_dir, "champion_model_forward_forecast.csv"))
ggplot2::ggsave(
  filename = file.path(output_dir, "historical_series.png"),
  plot = historical_plot,
  width = 12,
  height = 6,
  dpi = 300
)
ggplot2::ggsave(
  filename = file.path(output_dir, "backtest_forecast_comparison.png"),
  plot = backtest_plot,
  width = 12,
  height = 6,
  dpi = 300
)
writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "session_info.txt"))

cat(sprintf("\nAnalysis artifacts written to: %s\n", output_dir))