function future_prediction = Construct_Prediction(periods_count, coefficients_matrix, artificial_noise, hours_per_period, forecast_horizon, adjusted_data_series)
% Construct_Prediction - Constructs a prediction for future data using coefficients and noise
%   This function constructs a prediction for future data based on past data,
%   using a combination of coefficients and artificial noise.
%
%   Inputs:
%     - periods_count: Number of past periods used for prediction
%     - coefficients_matrix: Matrix containing coefficients obtained from past data
%     - artificial_noise: Artificial noise data added to improve prediction accuracy
%     - hours_per_period: Number of hours in each period of the data
%     - forecast_horizon: Number of future periods to forecast
%     - adjusted_data_series: Adjusted historical data series used in prediction
%
%   Output:
%     - future_prediction: Prediction for future data

    k = 1;
    % Adjust k value based on forecast horizon and periods_count
    if -periods_count + forecast_horizon * k >= periods_count
        k = k;
    else
        k = k + 1;
    end

    period_length = periods_count + 1;
    % Construct the matrix Q for prediction
    for i = 1:forecast_horizon * hours_per_period - periods_count
        if period_length == forecast_horizon + 1
            period_length = 1;
        end
        % Fill the Q matrix with appropriate coefficients
        for j = 1:periods_count + 1
            Q_Matrix(i, i + j - 1) = -coefficients_matrix(periods_count + 1 - j + 1, period_length);
        end
        period_length = period_length + 1;
    end
    % Slice the Q matrix to get the appropriate section
    Q_Matrix1 = Q_Matrix(:, periods_count + 1:end);
    % Calculate prediction using matrix inversion
    YY(periods_count + 1:(forecast_horizon * hours_per_period), 1) = inv(Q_Matrix1) * artificial_noise(periods_count + 1:(forecast_horizon * hours_per_period), 1);
    % Calculate predictions for the initial periods
    for i = 1:periods_count
        coefficients = flip(coefficients_matrix(2:periods_count + 1, i));
        YY(i, 1) = artificial_noise(i, 1) + dot(adjusted_data_series(k * forecast_horizon - periods_count + i:k * forecast_horizon - 1 + i), coefficients);
    end
    % Output the future prediction
    future_prediction = YY;
end
