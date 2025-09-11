library(tidyverse)
library(lubridate)
library(hms)
library(forecast)
library(tseries)
library(ggplot2)
library(broom)

#Load & preproces
raw <- read_csv(
  "london_smart_meter_with_timestamps.csv",
  col_types = cols(.default = col_double(), date = col_date())
)

#Remove last timestamp row
data <- raw[!is.na(raw$date), ]

#Reshape
data_long <- data |>
  pivot_longer(cols = -date, names_to = "time", values_to = "consumption") |>
  mutate(
    time = as_hms(time),
    datetime = as.POSIXct(date + time)
  ) |>
  drop_na()

#Daily average
daily_avg <- data_long |>
  group_by(date) |>
  summarise(daily_consumption = mean(consumption, na.rm = TRUE)) |>
  ungroup()

ts_daily <- ts(na.omit(daily_avg$daily_consumption), frequency = 7)

#Stationarity 
adf_result <- adf.test(ts_daily)
kpss_result <- kpss.test(ts_daily)

write_lines(c(
  paste("ADF Test p-value:", adf_result$p.value),
  paste("KPSS Test p-value:", kpss_result$p.value)
), "stationarity_tests.txt")

d_order <- ndiffs(ts_daily)
D_order <- nsdiffs(ts_daily)

ts_stationary <- ts_daily
if (d_order > 0) ts_stationary <- diff(ts_stationary, differences = d_order)
if (D_order > 0) ts_stationary <- diff(ts_stationary, lag = frequency(ts_daily), differences = D_order)

#ACF/PACF 
png("fig-acf-pacf.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
acf(ts_stationary, lag.max = 30, main = "ACF (stationary)")
pacf(ts_stationary, lag.max = 30, main = "PACF (stationary)")
dev.off()

#Spectral Estimation
spec_result <- spectrum(ts_daily, log = "yes", plot = FALSE)
dominant_freq <- spec_result$freq[which.max(spec_result$spec)]
dominant_period <- round(1 / dominant_freq)

cat("Dominant frequency:", dominant_freq, "\n")
cat("Dominant period (rounded):", dominant_period, "\n")
cat("Dominant period estimated from spectrum:", dominant_period, "\n")


png("fig-spectrum-log.png", width = 1100, height = 900)
spectrum(ts_daily, log = "yes", main = "Log-scaled Spectral Density")
abline(v = dominant_freq, col = "red", lty = 2)
text(dominant_freq, max(spec_result$spec), labels = sprintf("Period ≈ %.2f", 1 / dominant_freq), pos = 3, col = "blue")
dev.off()

write_csv(
  tibble(
    frequency = spec_result$freq,
    spectral_density = spec_result$spec,
    approx_period = 1 / spec_result$freq
  ),
  "spectral_density_log_scaled.csv"
)

#Manual Models
manual_models <- list(
  AR     = Arima(ts_daily, order = c(2, 0, 0)),
  MA     = Arima(ts_daily, order = c(0, 0, 2)),
  ARMA   = Arima(ts_daily, order = c(2, 0, 2)),
  ARIMA  = Arima(ts_daily, order = c(1, d_order, 1)),
  SARIMA = Arima(ts_daily, order = c(2, d_order, 2),
                 seasonal = list(order = c(2, D_order, 2), period = 7))
)

#SARIMA (Spectrum-Optimized)
#drawback if the dominantperiod not fit
if (is.finite(dominant_period) && dominant_period >= 2 && dominant_period <= 60) {
  seasonal_period <- round(dominant_period)
} else {
  seasonal_period <- 7
  warning("Invalid dominant period from spectral analysis. Using fallback period = 7.")
}

#fit SARIMA with spectrum
tryCatch({
  sarima_spec <- auto.arima(
    ts_daily,
    seasonal = TRUE,
    seasonal.period = seasonal_period,
    stepwise = FALSE,
    approximation = FALSE
  )
}, error = function(e) {
  message("Failed to fit spectrum-based SARIMA. Falling back to default seasonal model.")
  sarima_spec <<- Arima(
    ts_daily, 
    order = c(2, d_order, 2),
    seasonal = list(order = c(2, D_order, 2),
    period = seasonal_period))
})


#auto.arima model 
auto_model <- auto.arima(ts_daily, seasonal = TRUE)

#Time Series Regression
daily_avg <- daily_avg |>
  mutate(
    t = row_number(),
    weekday = factor(wday(date, label = TRUE, week_start = 1))
  )

lm_model <- lm(daily_consumption ~ t + weekday, data = daily_avg)

future_dates <- tibble(
  date = seq(max(daily_avg$date) + 1, by = "1 day", length.out = 30),
  t = (max(daily_avg$t) + 1):(max(daily_avg$t) + 30),
  weekday = factor(
    wday(seq(max(daily_avg$date) + 1, by = "1 day", length.out = 30),
         label = TRUE, week_start = 1),
    levels = levels(daily_avg$weekday))
)

pred_lm <- predict(lm_model, newdata = future_dates)

png("fig-forecast-tslm.png", width = 800, height = 400)
plot(ts(pred_lm, start = end(ts_daily)[1] + 1, frequency = 7),
     type = "l", col = "blue", lwd = 2,
     main = "Linear Regression Forecast (Trend + Weekday)",
     ylab = "kWh", xlab = "Time")
dev.off()

#Diagnostics
diagnostic_check <- function(model) {
  lb <- Box.test(residuals(model), lag = 20, type = "Ljung-Box")
  tibble(
    AIC = AIC(model),
    BIC = BIC(model),
    LB_pvalue = lb$p.value,
    RMSE = sqrt(mean(model$residuals^2, na.rm = TRUE))
  )
}

manual_diag <- map_dfr(manual_models, diagnostic_check, .id = "Model")
auto_diag <- diagnostic_check(auto_model) |> mutate(Model = "auto.arima()")
spec_diag <- diagnostic_check(sarima_spec) |> mutate(Model = "SARIMA_spectrum")
lm_rmse <- sqrt(mean(residuals(lm_model)^2))
lm_mae <- mean(abs(residuals(lm_model)))

model_diag <- bind_rows(manual_diag, auto_diag, spec_diag)
write_csv(model_diag, "model_diagnostics_updated.csv")

#Forecasting 
models_all <- c(
  manual_models,
  list(`auto.arima()` = auto_model, SARIMA_spectrum = sarima_spec)
)

fc_all <- map(models_all, forecast, h = 30)

walk2(fc_all, names(fc_all), \(fc, name) {
  ggsave(
    paste0("fig-forecast-", tolower(gsub("[.()]", "", name)), ".png"),
    autoplot(fc) +
      labs(title = paste(name, "30-day Forecast"), y = "kWh", x = "Time") +
      theme_minimal(),
    width = 8, height = 4
  )
})

#Evaluation 
eval_tbl <- tibble(
  Model = names(fc_all),
  RMSE  = map_dbl(fc_all, ~ sqrt(mean((.x$residuals)^2, na.rm = TRUE))),
  MAE   = map_dbl(fc_all, ~ mean(abs(.x$residuals), na.rm = TRUE)),
  AIC   = map_dbl(models_all, AIC),
  BIC   = map_dbl(models_all, BIC)
)

eval_tbl <- bind_rows(
  eval_tbl,
  tibble(Model = "TS_Regression", RMSE = lm_rmse, MAE = lm_mae,
         AIC = AIC(lm_model), BIC = BIC(lm_model))
)

write_csv(eval_tbl, "forecast_model_comparison_final.csv")

#Residual Checks
for (name in names(models_all)) {
  png(paste0("fig-resid-", tolower(name), ".png"), width = 700, height = 400)
  checkresiduals(models_all[[name]], main = paste(name, "Residual Check"))
  dev.off()
}
