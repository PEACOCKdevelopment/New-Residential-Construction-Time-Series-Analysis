# New Residential Construction Time Series Analysis

A reproducible, audit-friendly forecasting project that examines monthly U.S. new residential construction starts and turns raw observations into decision-ready forecast artifacts.

## Project Overview

This repository analyzes monthly housing-start data (Aug 1985-Feb 2025) using a compact production workflow in R and R Markdown. The project focuses on practical model evaluation, not only in-sample fit.

Key objectives:
- Understand long-term cyclical and seasonal behavior in residential construction activity.
- Benchmark additive and multiplicative Holt-Winters variants against baseline models.
- Validate real forecasting performance with a strict out-of-sample backtest.
- Export transparent artifacts for reporting, review, and reproducibility.

## What Makes It Interesting

Housing starts are a high-signal macro indicator, but seasonal structure can easily mislead model selection. This project highlights a common real-world finding: the model with the best training fit is not always the best forecaster. By separating train and test performance, it demonstrates why robust backtesting matters in economic time-series work.

## Methodology At A Glance

- Data quality controls: date parsing checks, NA checks, duplicate-period detection, continuity checks, and positivity validation.
- Seasonal-trend decomposition with STL for structural interpretation.
- Model set:
	- Holt-Winters Additive (Auto)
	- Holt-Winters Additive (Fixed)
	- Holt-Winters Multiplicative (Auto)
	- Holt-Winters Multiplicative (Fixed)
	- Naive baseline
	- Seasonal Naive baseline
- Evaluation metrics: MAE, MSE, RMSE, MAPE, and AMAPE (sMAPE label in script output).
- Residual diagnostics with Ljung-Box testing.
- Champion-model selection based on out-of-sample test MAE.

## Analysis Results

- Sample size: 475 monthly observations (Aug 1985-Feb 2025).
- Best in-sample model: HW Multiplicative (Auto), MAE 7.24, RMSE 9.38, MAPE 7.09%.
- Best out-of-sample model (12-month test): HW Additive (Fixed), MAE 6.46, RMSE 7.69, MAPE 5.92%, AMAPE 5.69%.
- Benchmark comparison on test window:
	- HW Additive (Fixed): MAE 6.46
	- Naive: MAE 7.23
	- Seasonal Naive: MAE 7.50
- Test MAE improvement of the champion model:
	- 10.75% better than Naive
	- 13.92% better than Seasonal Naive
- Residual diagnostics: Ljung-Box p-values are low across models, suggesting remaining autocorrelation not fully captured by Holt-Winters specifications.
- 12-month forward forecast (champion model): point forecasts range from 97.02 to 124.23 thousand units; first month (Mar 2025) forecast is 112.32 and last month (Feb 2026) forecast is 105.96.

## Outputs

Generated artifacts are saved in `analysis_outputs/`, including:
- `model_metrics.csv`
- `holt_winters_model_metrics.csv`
- `model_residual_diagnostics.csv`
- `champion_model_forward_forecast.csv`
- `session_info.txt`

The repository also contains an HTML report with narrative interpretation and visual diagnostics.

## Live Report (RPubs)

Read the published report here:

https://rpubs.com/timothypawelczyk/new-construction-analysis

## Reproducibility

Run the production script:

```bash
Rscript construction_data_analysis.R
```

Render the report:

```bash
Rscript -e "rmarkdown::render('construction_data_analysis.Rmd', output_file = 'construction_data_analysis.html')"
```

## Repository Structure

- `construction_data_analysis.R` - production analysis pipeline
- `construction_data_analysis.Rmd` - narrative report source
- `construction_data_analysis.html` - rendered report
- `construction_data.csv` - input dataset
- `analysis_outputs/` - exported metrics and forecast artifacts

## Author

Tymoteusz Pawelczyk
