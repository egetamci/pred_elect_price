% Clear all previous variables and data
clear all;

% Load the CSV file containing price data
opts = detectImportOptions('... /prices.csv');
opts.VariableNamingRule = 'preserve';
data_table = readtable('... /prices', opts);

% Extract and clean price data
prices = data_table.PTF;
cleaned_prices = str2double(strrep(strrep(data_table.PTF, '.', ''), ',', ''));
initial_data  = cleaned_prices / 100000;
initial_data = initial_data(1:end-120,1);
initial_data_length = length(initial_data);

% Define the forecast horizon
forecast_horizon = 24;

% Detrend the data
detrend_data = detrend(initial_data);
trend = initial_data - detrend_data;

% Extract an initial window for further processing
initial_window = detrend_data(1:initial_data_length - forecast_horizon, 1);
initial_window_length = length(initial_window);

% Calculate initial and prediction hours
initial_hours = initial_window_length / forecast_horizon;
prediction_hours = initial_data_length / forecast_horizon;

% Generate sine-cosine data
time_points = 1:initial_window_length;
sine_cosine_data = Sinusoidal(initial_window, initial_window_length, forecast_horizon);
adjusted_data = initial_window - sine_cosine_data';

% Calculate auto-period and period averages
auto_period = Covariance(adjusted_data, 100);
for i = 1:forecast_horizon
    period_sum = sum(adjusted_data(i:forecast_horizon:forecast_horizon*(initial_hours-1)+i,1));
    period_average(i, 1) = period_sum / initial_hours;
end

% Adjust data with period averages
for i = 1:forecast_horizon
    for j = 1:initial_hours
        adjusted_data_periodic(i+(j-1)*forecast_horizon, 1) = adjusted_data(i+(j-1)*forecast_horizon, 1) - period_average(i);
    end
end

% Copy adjusted data and period averages
adjusted_data_copy = adjusted_data_periodic;
period_copy = period_average;

% Perform Akaike Information Criterion (AIC) calculation
[Aic_qp1, num_periods] = AIC_Parameter_Selection(forecast_horizon, adjusted_data, initial_hours);
aa_matrix = -ones(num_periods+1, forecast_horizon);
coefficients = ForecastCoeffs(num_periods, adjusted_data_copy, initial_hours, forecast_horizon);
aa_matrix(2:end, :) = coefficients';
noise = APP_noise(num_periods, aa_matrix, adjusted_data_copy, initial_hours, forecast_horizon);

% Perform iterations for noise generation and prediction construction
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
predicted_values = forecasting(forecast_horizon, period_copy, sine_cosine_data', initial_hours, prediction_hours, final_data);
predicted_values = predicted_values + trend_difference;
initial_data = detrend_data + trend_difference;

% Calculate mean prediction error
mean_prediction_error = 0;
for i = 1:forecast_horizon
    prediction_error(i) = (abs(initial_data(end-i+1, 1) - predicted_values(end-i+1, 1)))^2;
    mean_prediction_error = mean_prediction_error + prediction_error(i);
end
mean_prediction_error = (mean_prediction_error)/ forecast_horizon
% Plot the results
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

