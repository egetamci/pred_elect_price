function [AIC_values, optimal_pp] = AIC_Parameter_Selection(parameters, input_data, data_length)
% Parameter_Selection_AIC - Selects parameters using Akaike Information Criterion (AIC) for α-stable PAR model
%   [AIC_values, optimal_pp] = Parameter_Selection_AIC(parameters, input_data, data_length)
%   This function performs parameter selection for an α-stable PAR model
%   using the Akaike Information Criterion (AIC).
%
%   Inputs:
%     - parameters: Number of parameters
%     - input_data: Input data
%     - data_length: Length of the data
%
%   Outputs:
%     - AIC_values: Vector containing AIC values for different parameter selections
%     - optimal_pp: The value of pp corresponding to the minimum AIC value
    
    % Initialize minimum AIC value
    min_AIC = 10^(10);
    % Initialize variable to store the value of pp corresponding to minimum AIC
    optimal_pp_index = 0;
    % Set the maximum value of pp
    max_pp = 15;
    
    % Loop through different values of pp
    for pp = 1:max_pp
        % Initialize coefficient matrix
        coefficients_matrix = -ones(pp + 1, parameters);
        % Calculate coefficients using Calculate_Coefficients function
        coefficients = Calculate_Coefficients(pp, input_data, data_length, parameters);
        % Update the coefficient matrix
        coefficients_matrix(2:end, :) = coefficients';
        % Calculate noise using Calculate_Noise function
        noise = Calculate_Noise(pp, coefficients_matrix, input_data, data_length, parameters);
        % Estimate parameters using stblfit function
        alpha = stblfit(noise);
        % Calculate AIC using Calculate_AIC function
        AIC_values(pp, 1) = Calculate_AIC(noise, alpha(1, 1), 0, 1, 0, parameters, pp);
        
        % Update minimum AIC value and corresponding pp value if necessary
        if (AIC_values(pp, 1) < min_AIC)
            min_AIC = AIC_values(pp, 1);
            optimal_pp_index = pp;
        end
    end
    % Store the value of pp corresponding to minimum AIC
    optimal_pp = optimal_pp_index;
end

function AIC = Calculate_AIC(noise, alpha, beta, gamma, Mu, period, pp)
    % Calculate the probability density function of the α-stable distribution
    alpha_stable_pdf = stblpdf(noise, alpha, beta, gamma, Mu);
    % Take the natural logarithm of the probability density function
    alpha_stable_log_pdf = log(alpha_stable_pdf);
    % Calculate the AIC value
    AIC = -2 * (sum(alpha_stable_log_pdf)) + 2 * period * (pp + 1);
end

function coefficients = Calculate_Coefficients(pp, data, data_length, parameters)
    % Initialize matrices for covariance and R values
    covariance_matrix = zeros(pp, pp);
    R_vector = zeros(pp, 1);
    
    % Loop through iterations
    for iter = 1:parameters
        % Compute covariance matrix using Compute_Covariance function
        covariance_matrix(:,:) = Compute_Covariance(iter, pp, data, data_length, parameters);
        % Compute R vector using Compute_R function
        R_vector(:, 1) = Compute_R(pp, data, data_length, parameters, iter);
        % Compute coefficients using matrix inversion
        coefficients(iter, :) = inv(covariance_matrix) * R_vector;
    end
end

function covariance_matrix = Compute_Covariance(lag_index, lag_count, signal_data, data_length, parameter_count)
% Compute_Covariance - Computes the covariance matrix for a given signal.
%   This function computes the covariance matrix for a given signal up to a
%   specified lag count.
%
%   Inputs:
%     - lag_index: Index of the lag to compute covariance for
%     - lag_count: Total number of lags for covariance computation
%     - signal_data: Input signal data
%     - data_length: Length of the input data
%     - parameter_count: Number of parameters
%
%   Output:
%     - covariance_matrix: Covariance matrix

    % Initialize covariance matrix
    covariance_matrix = zeros(lag_count, lag_count);
    
    % Loop through each element of the covariance matrix
    for row = 1:lag_count
        for col = 1:lag_count
            % Compute each element using Lambda function
            covariance_matrix(row, col) = Lambda(signal_data, data_length, lag_index - col, row - col, parameter_count);
        end
    end
end

