
Electricity Price Forecasting using Generalized Yule-Walker Method
This repository contains MATLAB code implementing the Generalized Yule-Walker Method for forecasting electricity prices in the Turkish market. Basically, It is used that the PAR series outperforms with SαS which is the special class of α stable distributions, as an expansion of the well-defined AR series.

Main Steps:

1] Data Loading and Preprocessing:
* Load the CSV file containing price data.
* Clean and preprocess the data.
  
2] Detrending and Initial Window Extraction:
* Detrend the data to remove any trend components.
* Extract an initial window for further processing.

3] Generating Sine-Cosine Data:
* Generate sine-cosine data to capture periodic patterns in the data.

4] Adjusting Data with Period Averages:
* Calculate auto-period and period averages.
* Adjust data with period averages to remove periodic components.
  
5] AIC Calculation:
* Perform Akaike Information Criterion (AIC) calculation to determine the optimal model order.
  
6] Generating Noise for Prediction:
* Generate noise for prediction construction using stable distribution parameters.
  
7] Constructing Predictions:
* Use coefficients and noise to construct predictions for future data points.
  
8] Combining Predictions with Trend:
* Combine predictions with the trend component to generate the final forecast.
  
9] Forecast Evaluation:
 * Evaluate the accuracy of the forecast using appropriate metrics.

References

F. A. Savacı & E. Tamcı, 2023, ”Forecasting Turkish Electricity Prices with Alpha-Stable Periodic Autoregressive Time Series”, ELECO 2023 14th Int. Conf. Electrical and Electronics Engineering, 1-5

