

Forecasting Turkey Electricity Prices by Alpha-Stable Periodic AutoRegressive Time Series 
			

In this study, Electricity Price Forecasting is conducted using the Generalized Yule-Walker Method, with MATLAB code available in the repository for implementation in the Turkish market [1]. Our findings reveal that the periodic autoregressive (PAR) time series outperforms the SS, a special class of $\alpha$-stable distributions, extending the well-established autoregressive (AR) series. The choice of the "Generalized Yule-Walker Method" is motivated by the absence of finite variance in random variables within the classic Yule-Walker model, which poses challenges in coefficient evaluation for PAR models with i.i.d. random variables from -stable distribution [2].

Table of Contents

Model
Forecasting Process:
Data Loading and Preprocessing:
 Detrending and Initial Window Extraction:
 Generating Sine-Cosine Data:
Adjust data with period averages to remove periodic components.
Evaluating the Order (p) of the Model by Akaike Information Criterion (AIC):
Perform iterations for noise generation and prediction construction
Ploting the Forecasted Electricity Prices:


     	                                                   (1)
