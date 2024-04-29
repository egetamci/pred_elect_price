function forecast_signal = forecasting(p, period_array, seasonal_component, total_data_length, forecast_length, observed_signal)
    % Compute the starting index for the seasonal component
    start_index = total_data_length - forecast_length;
    
    % Initialize forecast signal array
    forecast_signal = zeros(p * total_data_length, 1);
    
    % Generate forecast signal
    for i = 1:p
        for j = 1:total_data_length
            forecast_signal(i + (j - 1) * p, 1) = observed_signal(i + (j - 1) * p, 1) + period_array(i, 1);
        end
    end
    
    % Adjust forecast signal with seasonal component
    forecast_signal(1:p * start_index, 1) = forecast_signal(1:p * start_index, 1) + seasonal_component;
    
    % Compute the average seasonal component
    average_seasonal_component = zeros(p, 1);
    for i = 1:forecast_length
        average_seasonal_component = average_seasonal_component + seasonal_component(1 + (p * (i - 1)):p * i, 1);
    end
    average_seasonal_component = average_seasonal_component ./ forecast_length;
    
    % Duplicate the average seasonal component for missing data
    adjusted_seasonal_component = zeros(p, 1);
    for k = 1:start_index
        adjusted_seasonal_component(1 + (k - 1) * p:k * p, 1) = average_seasonal_component;
    end
    
    % Add adjusted seasonal component to the forecast signal
    forecast_signal(p * start_index + 1:p * total_data_length, 1) = forecast_signal(p * start_index + 1:p * total_data_length, 1) + adjusted_seasonal_component;
end
