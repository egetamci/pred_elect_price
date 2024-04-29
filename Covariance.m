function covariance_sequence = Covariance(X, num_lags)
% Covariance - Computes the covariance sequence of a given signal.
%   This function computes the covariance sequence of a given signal vector
%   up to a specified number of lags.
%
%   Inputs:
%     - X: Input signal vector
%     - num_lags: Number of lags for computing covariance
%
%   Output:
%     - covariance_sequence: Covariance sequence

    signal_length = length(X);
    total_signal_abs_sum = sum(abs(X));
    for lag = 1:num_lags
        lag_sum = 0;
        for k = lag + 1:signal_length
            lag_product = X(k) * sign(X(k - lag));
            lag_sum = lag_sum + lag_product;
        end
        covariance_sequence(lag) = lag_sum;
    end
    covariance_sequence = covariance_sequence ./ total_signal_abs_sum;
end

