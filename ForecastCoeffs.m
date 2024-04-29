
function coeff5 = ForecastCoeffs(num_periods_2, adjusted_data_copy_2, initial_hours_2, forecast_horizon_2)
% Calculates coefficients for forecasting based on input data and parameters.
% Inputs:
%   - num_periods_2: Number of periods
%   - adjusted_data_copy_2: Adjusted data copy
%   - initial_hours_2: Initial hours
%   - forecast_horizon_2: Forecast horizon
% Output:
%   - coeff5: Coefficients matrix for forecasting
Lg_m1c=zeros(num_periods_2,num_periods_2);
Rg_m1c=zeros(num_periods_2,1);
for ii = 1:forecast_horizon_2
Lg_m1c = calculate_covariance(ii, num_periods_2, adjusted_data_copy_2, initial_hours_2,forecast_horizon_2);
Rg_m1c = RM(num_periods_2, adjusted_data_copy_2, initial_hours_2,forecast_horizon_2, ii);
coeff5(ii, :) = inv(Lg_m1c).*Rg_m1c;
end
end
function RM1=R_M(pp22,X22,N22,T22,vv2)
for j2=1:pp22
RM1(j2)=Lamda_M(X22,T22,vv2,j2,N22);
end
end

function covariance_matrix = calculate_covariance(ss1, num_periods_2, adjusted_data_copy_2, initial_hours_2, forecast_horizon_2)
% Calculates a covariance-like matrix based on the provided parameters.
% Inputs:
%   - ss1: Some parameter
%   - num_periods_2: Number of periods
%   - adjusted_data_copy_2: Adjusted data copy
%   - initial_hours_2: Initial hours
%   - forecast_horizon_2: Forecast horizon
% Output:
%   - covariance_matrix: Computed covariance-like matrix

    for i3 = 1:num_periods_2
        for j3 = 1:num_periods_2
            % Calculate each element of the covariance-like matrix
            covariance_matrix(i3, j3) = calculate_Lamda_coeff(adjusted_data_copy_2, forecast_horizon_2, ss1 - j3, i3 - j3, initial_hours_2);
        end
    end
end

function RM1 = RM(num_periods_2, adjusted_data_copy_2, initial_hours_2, forecast_horizon_2, ii_2)
%   RM function estimates coefficients for an α-stable PAR model and returns
%   the resulting coefficients vector.
%   Inputs:
%     - num_periods_2: Number of periods
%     - adjusted_data_copy_2: Adjusted data series
%     - initial_hours_2: Initial hours
%     - forecast_horizon_2: Forecast horizon
%     - ii_2: Parameter for estimation
%
%   Output:
%     - RM1: Resulting coefficients vector 
% Initialize the resulting coefficients vector
    RM1 = zeros(num_periods_2, 1);
    
    % Loop through each period and estimate coefficients using Lamda_M function
    for j_2 = 1:num_periods_2
        RM1(j_2) = Lamda_M(adjusted_data_copy_2, forecast_horizon_2, ii_2, j_2, initial_hours_2);
    end
end


function lamda = Lamda_M(adjusted_data_copy_2, forecast_horizon_2, v11, k11, initial_hours_2)
% Lamda_M - Estimate Parameters for α-stable PAR Model
%   Lamda_M function estimates the parameters for an α-stable PAR model and 
%   returns the coefficients of the resulting time series.
%
%   Inputs:
%     - adjusted_data_copy_2: Adjusted data series
%     - forecast_horizon_2: Forecast horizon
%     - v11, k11: Parameters for estimation
%     - initial_hours_2: Initial hours
%
%   Output:
%     - lamda: Coefficients of the resulting time series
%
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
    
    % Compute the coefficient
    lamda = ab_u / ab_d;
end

