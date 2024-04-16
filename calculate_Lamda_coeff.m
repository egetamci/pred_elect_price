function lamda = calculate_Lamda_coeff(adjusted_data_copy_2, forecast_horizon_2, v11, k11, initial_hours_2)
    % Adjust parameters if k11 is negative
    if k11 < 0
        v11 = v11 - k11;
        k11 = -k11;
    end
    
    % Compute r11 based on parameters
    if (v11 > 0) || (v11 - k11 > 0)
        r11 = initial_hours_2 - 1;
    else
        r11 = initial_hours_2;
    end
    
    % Compute l11 based on parameters
    if (v11 <= 0) || (v11 - k11 <= 0)
        l11 = 1;
    else
        l11 = 0;
    end
    
    % Compute ab and ab1 based on adjusted data
    ab = adjusted_data_copy_2(l11 * forecast_horizon_2 + v11 : forecast_horizon_2 : r11 * forecast_horizon_2 + v11);
    ab1 = sign(adjusted_data_copy_2(l11 * forecast_horizon_2 + v11 - k11 : forecast_horizon_2 : r11 * forecast_horizon_2 + v11 - k11));
    
    % Compute the numerator and denominator for the coefficient estimation
    ab_u = dot(ab, ab1);
    ab_d = sum(abs(adjusted_data_copy_2(l11 * forecast_horizon_2 + v11 - k11 : forecast_horizon_2 : r11 * forecast_horizon_2 + v11 - k11)));
    
    % Compute the coefficient for the covariance matrix element
    lamda = ab_u / ab_d;
end
