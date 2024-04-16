function YY1 = Construct_Prediction(num_periods_1, aa_matrix_1, artificial_noise_1, hours_1, forecast_horizon_1, adjusted_data_copy_1)
  % Construct_Prediction - Construct Prediction using Coefficients and Noise
%   Construct_Prediction function constructs the prediction for the future data
%   using coefficients and noise.
%
%   Inputs:
%     - num_periods_1: Number of periods
%     - aa_matrix_1: Coefficients matrix
%     - artificial_noise_1: Artificial noise data
%     - hours_1: Hours
%     - forecast_horizon_1: Forecast horizon
%     - adjusted_data_copy_1: Adjusted data series
%
%   Output:
%     - YY1: Constructed prediction for the future data
%
  k_pp = 1;
    if -num_periods_1 + forecast_horizon_1 * k_pp >= num_periods_1
        k_pp = k_pp;
    else
        k_pp = k_pp + 1;
    end
    
    peri = num_periods_1 + 1;
    for i5 = 1:forecast_horizon_1 * hours_1 - num_periods_1
        if peri == forecast_horizon_1 + 1
            peri = 1;
        end
        for j3 = 1:num_periods_1 + 1
            Q_Mat(i5, i5 + j3 - 1) = -aa_matrix_1(num_periods_1 + 1 - j3 + 1, peri);
        end
        peri = peri + 1;
    end
    Q_Mat1 = Q_Mat(:, num_periods_1 + 1:end);
    YY(num_periods_1 + 1:(forecast_horizon_1 * hours_1), 1) = inv(Q_Mat1)* artificial_noise_1(num_periods_1 + 1:(forecast_horizon_1 * hours_1),1);
    for i22 = 1:num_periods_1
        aa_f2 = flip(aa_matrix_1(2:num_periods_1 + 1, i22));
        YY(i22, 1) = artificial_noise_1(i22, 1) + dot(adjusted_data_copy_1(k_pp * forecast_horizon_1- num_periods_1 + i22:k_pp * forecast_horizon_1 - 1 + i22), aa_f2);
    end
    YY1 = YY;
end
