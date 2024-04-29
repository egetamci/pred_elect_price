 # Forecasting Turkey Electricity Prices by Alpha-Stable Periodic AutoRegressive Time Series 
			
In this study, Electricity Price Forecasting is conducted using the Generalized Yule-Walker Method, with MATLAB code available in the repository for implementation in the Turkish market [1]. Our findings reveal that the periodic autoregressive (PAR) time series outperforms the $S\alpha S$, a special class of $\alpha$-stable distributions, extending the well-established autoregressive (AR) series. The choice of the "Generalized Yule-Walker Method" is motivated by the absence of finite variance in random variables within the classic Yule-Walker model, which poses challenges in coefficient evaluation for PAR models with i.i.d. random variables from $\alpha$-stable distribution [2].

## Table of Contents
* Model
* Forecasting Process:
* Data Loading and Preprocessing:
* Detrending and Initial Window Extraction:
* Generating Sine-Cosine Data:
* Adjust data with period averages to remove periodic components.
* Evaluating the Order (p) of the Model by Akaike Information Criterion (AIC)
* Perform iterations for noise generation and prediction construction
* Ploting the Forecasted Electricity Prices:

**Model**


   $$X_t-\Phi_1(t) X_{t-1}-\cdots-\Phi_p(t) X_{t-p}=\xi_t$$

Here,  $\Phi_i(t) \in \mathbb{R}, ~i = 1, \ldots, p,$  are periodic in  with the same period $T$ and the length of the sequence $\{X_t\}$  is equal to .

**Forecasting Proces**

* Data Loading and Preprocessing:

```

clear all;

opts = detectImportOptions('/Users/egetamci/Downloads/prices.csv')

opts.VariableNamingRule = 'preserve';

data_table = readtable('/Users/egetamci/Downloads/prices', opts);
```

*   Clean and preprocess the data.
   
```   
prices = data_table.PTF;

cleaned_prices = str2double(strrep(strrep(data_table.PTF, '.', ''), ',', ''));

initial_data = cleaned_prices / 100000;

initial_data = initial_data(1:end-120,1);

initial_data_length = length(initial_data);
 ```
**Detrending and Initial Window Extraction:**

* Detrend the data to remove any trend components.
```
detrend_data = detrend(initial_data);

trend = initial_data - detrend_data;
```
* Extract an initial window for further processing.
```
initial_window = detrend_data(1:initial_data_length - forecast_horizon, 1);

initial_window_length = length(initial_window);
```
**Generating Sine-Cosine Data:**

* Generate sine-cosine data to capture periodic patterns in the data.
```
time_points = 1:initial_window_length;

sine_cosine_data = Sinusoidal(initial_window, initial_window_length, forecast_horizon);

adjusted_data = initial_window - sine_cosine_data';
```

Then to estimate the period  by a meaure which is used for independence is  function of the stationary process  $X_t$ for lag $k$ is defined as the covariation of random variables $X_{t}$ and $X_{t-k}$ , which is $CV(X_t, X_{t-k})$.

Estimation is defined over normalized autocovariation function . Namely, it is propose to estimate the value:

The normalized autocovariation for process $X_{t}$ for lag $k$ is given by:

$$
 NCV(X_t,X_{t-k})=\frac{E(X_t sign(X_{t-k}))}{E|X_{t-k}|}.
$$

Where $NCV(X_t,X_{t-k})=\frac{CV(X_t , X_{t-k})}{{\sigma}^{\alpha}{X{t-k}}}$ and $\sigma_{X_{t-k}}$ is the scale paramere of the $X_{t-k}$  random variable.
Estimator of the NCV for a single trajectory a sample   being a realization of a stationary process  takes the form:

$$
\hat{NCV}(X_t, X_{t-k})=\frac{ {\sum}^r_{t=l} x_t\text{sign}(X_{t-k})}{{\sum}^N_{t=l}\|x_t\|}
$$

$$
l= \max(1, 1 + k)
$$
$$
r= \min(N, N + k).
$$
```
auto_period = Covariance(adjusted_data, 100);
```


Based on the graphic provided above evaluated by Covariance function, 24 hours is estimated as the period $T.$

Then, the periodic mean $\mu(v)$  is estimated by the period  $T$  :

$$
\hat{\mu}(v)=\frac{1}{N}\sum^{N-1}_{n=0} X_nT+v, \quad v=1,\ldots,T
$$




where $N =\frac{L}{T}$

Then, the mean is distracted the mean from the data to use the method.

$$
Y_{nT+v}=X_{nT+v}-\hat{\mu}(v)
$$

where $Y_{n}$ is the zero mean series.

**Adjust data with period averages to remove periodic components**

* Calculate auto-period and period averages.
```
auto_period = Covariance(adjusted_data, 100);

 for i = 1:forecast_horizon

     period_sum = sum(adjusted_data(i:forecast_horizon:forecast_horizon*(initial_hours-1)+i,1));

     period_average(i, 1) = period_sum / initial_hours;

 end
```
 * Adjust data with period averages to remove periodic components.
```
   for i = 1:forecast_horizon

     for j = 1:initial_hours

        adjusted_data_periodic(i+(j-1)*forecast_horizon, 1) = adjusted_data(i+(j-1)*forecast_horizon, 1) - 
        period_average(i);

     end

   end
```

Then, given AR model can be represented by considering $Y_{n}$ variable as follows:

$$
Y_{nT+v} - \sum_{i=1}^{p} \phi_i(nT+v) Y_{nT+v-i} = \xi_{nT+v}
$$

**Evaluating the Order (p) of the Model by Akaike Information Criterion (AIC):**

Another necessary issue is fitting the data to the model efficiency. To find out the order in the Model, the AIC Akaike information criterion (AIC) is employed for $\alpha$-stable process outlined below.

* Perform Akaike Information Criterion (AIC) calculation to determine the optimal model order.

$$
AIC(p) = -2 \sum_{i=1 }^{NT} \log(p_{S_{\alpha,\sigma,\beta,\mu}}(\xi_i)) + 2T(p+1).
$$

![AIC](https://github.com/egetamci/pred_elect_price/assets/160476027/c6fd4ffc-1759-4e8a-b469-008b5f3be9ed)



Here $p_{S_{\alpha,\sigma,\beta,\mu}}$ is the density function of the stable distribution. The smallest value of the AIC(p) function gives the most appropriate order for the model. The most optimal order is given for the smallest value of AIC(p).

```
[Aic_qp1, num_periods] = Akai_Inf_Cri(forecast_horizon, adjusted_data, initial_hours);
```
Other stages resume as known as in AR case, we multiply the above equation by vector $[S_{nT+v-1},...,S_{nT+v-p}]$  wherever $S_k = sign(X_k)$ , and obtain the system of p equations. Within the next step, the expectation is calculated either and take into account that the autocovariation is periodic with period $T$ :

$$
E [Y_v S_{v-1} - \sum^p_{i=1}\Phi_i(v) \mathbb{E}[X_{v-i}S_{v-1}]] = E[\xi_{nT+v}]
$$

.

.

.

$$
E [Y_v S_{v-p} - \sum^p_{i=1}\Phi_i(v) \mathbb{E}[X_{v-i}S_{v-p}]] = E[\xi_{nT+v}]
$$

To obtain the autocovariation,system is divided by $[E|X_{v-1}|, . . . E|X_{v-p}|]^{'}$  and takes
following form:

$$
\lambda_v^{(1)} - \sum^p_{i=1}\Phi_i(v)\lambda_v^{(1-i)} = 0
$$

.

.

.

$$
\lambda_v^{(p)} - \sum^p_{i=1}\Phi_i(v)\lambda_{v-i}^{(p-i)} = 0
$$


where $\lambda_v(k) = \frac{E[Y_v S_{v-k}]}{E[|X_{v-k}|]}, \quad v = 1, \ldots, T$ and $k = 1,...,p$.
 
Conclusively,  the parameters are rearranged to be the solution for the matrix system below:

$$
\lambda_v = \Lambda_v \Phi_v
$$

where $\lambda_v$  and $\Phi_v$  are vectors of length p defined as:

$$
\lambda_v =[\lambda_v(1)~...~  \lambda_v(p)]^T
$$

$$
\Phi_v =[\Phi_1(v)~...~ \Phi_p(v)]^T.
$$

The matrix $\Lambda_v$ has size $p\times p$ and its $(i,j)$ th part is given by:

$$
\Lambda_v(i,j)=\lambda_{v-i}(i-j)
$$

Here, the definition of the measures depends on the type of the data. Estimation of the autocovariance function is remodeled according to period. Besides, the estimator of auto variation is evaluated as follows wherever $v = 1, . . . , T.$ 

Finally, necessary parameters within the $\alpha$-stable AR model  can be calculated via $\hat{\Phi}^{-1}=\hat{\Lambda}^{-1}_v\hat{\lambda}^{-1}_v$ for the construction of the Model (1).

```
aa_matrix = -ones(num_periods+1, forecast_horizon);

coefficients = ForecastCoeffs(num_periods, adjusted_data_copy,initial_hours,forecast_horizon);

aa_matrix(2:end, :) = coefficients';
```
To ensure that the residuals are independent, a graphical test based on empirical autocovariation is suggested


```
noise = Noise_extraction(num_periods, aa_matrix, adjusted_data_copy, initial_hours, forecast_horizon);
auto_period_noise = Covariance(noise, 100);
```

**Perform iterations for noise generation and prediction construction**

```
auto_period_noise = Covariance(noise, 100);
coefficients_copy = coefficients;
maximum_noise = max(abs(noise));
stable_parameters = stblfit(noise);
artificial_noise = zeros((prediction_hours-initial_hours)*forecast_horizon, 1);
prediction_iterations = 1000;
forecast_data = zeros((prediction_hours-initial_hours)*forecast_horizon, 1);

for i = 1:prediction_iterations
    generated_noise = stblrnd(stable_parameters(1), 0, 1, 0, 1000, 1);
    artificial_noise(1:(prediction_hours-initial_hours)*forecast_horizon,1) = generated_noise(2*forecast_horizon+1:2*forecast_horizon+(prediction_hours-initial_hours)*forecast_horizon,1);
    noise_magnitude = max(abs(artificial_noise));
    artificial_noise = (maximum_noise / noise_magnitude) .* artificial_noise;
    prediction_array = Construct_Prediction(num_periods, aa_matrix, artificial_noise, prediction_hours-initial_hours, forecast_horizon, adjusted_data_copy);
    forecast_data = forecast_data + prediction_array;
end

% Average the forecast data over prediction iterations
forecast_data = forecast_data / prediction_iterations;

% Predict the linear trend of the data
p = polyfit(1:initial_data_length-forecast_horizon, detrend_data(1:end-forecast_horizon,1), 1);
x = 1:initial_data_length;
trend_difference = p(1).*x+p(2);

% Combine forecast data with trend
final_data(1:(initial_hours)*forecast_horizon, 1) = adjusted_data_periodic;
final_data(initial_hours*forecast_horizon+1:prediction_hours*forecast_horizon) = forecast_data;
predicted_values = forecasting_t1(forecast_horizon, period_copy, sine_cosine_data', initial_hours, prediction_hours, final_data);
predicted_values = predicted_values + trend_difference;
initial_data = detrend_data + trend_difference;

% Calculate mean prediction error
mean_prediction_error = 0;
for i = 1:forecast_horizon
    prediction_error(i) = (abs(initial_data(end-i+1, 1) - predicted_values(end-i+1, 1)))^2;
    mean_prediction_error = mean_prediction_error + prediction_error(i);
end
mean_prediction_error = (mean_prediction_error)/ forecast_horizon
```
**Ploting the Forecasted Electricity Prices:**

Finally, the one-day ahead prediction result of electricity prices is obtained
```
figure;
plot(initial_data(end-24:end, 1), 'g', 'LineWidth', 2); % Original data
hold on;
plot(predicted_values(end-24:end, 1), 'r', 'LineWidth', 2); % Predicted values

% Add labels and title
xlabel('Time');
ylabel('Value');
title('Comparison of Original Data and Predicted Values');

% Add legend
legend('Original Data', 'Predicted Values', 'Location', 'northwest');
```

**References:**

1 ] F. A. Savacı & E. Tamcı, 2023, ”Forecasting Turkish Electricity Prices with Alpha-Stable Periodic Autoregressive Time Series”, ELECO 2023 14th Int. Conf. Electrical and Electronics Engineering, 1-5.

2 ] T. M. G. J. Kruczek P, Wylomanska A, “The modified yulewalker method for a-stable time series models,” Physica A:Stat Mech Appl, vol. 469, p. 588–603, 2017.






