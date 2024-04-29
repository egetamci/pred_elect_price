function lamda_coefficient = calculate_Lamda_coefficient(adjusted_data, forecast_horizon, v, k, initial_hours)
    % Adjust parameters if k is negative
    if k < 0
        v = v - k;
        k = -k;
    end
    
    % Compute r based on parameters
    if (v > 0) || (v - k > 0)
        r = initial_hours - 1;
    else
        r = initial_hours;
    end
    
    % Compute l based on parameters
    if (v <= 0) || (v - k <= 0)
        l = 1;
    else
        l = 0;
    end
    
    % Extract subsequence from adjusted data based on parameters
    subsequence = adjusted_data(l * forecast_horizon + v : forecast_horizon : r * forecast_horizon + v);
    
    % Extract signs subsequence from adjusted data based on parameters
    signs_subsequence = sign(adjusted_data(l * forecast_horizon + v - k : forecast_horizon : r * forecast_horizon + v - k));
    
    % Compute numerator and denominator for coefficient estimation
    numerator = dot(subsequence, signs_subsequence);
    denominator = sum(abs(adjusted_data(l * forecast_horizon + v - k : forecast_horizon : r * forecast_horizon + v - k)));
    
    % Compute the coefficient for the covariance matrix element
    lamda_coefficient = numerator / denominator;
end
