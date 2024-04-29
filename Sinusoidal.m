function sinusoidal_signal = Sinusoidal(input_array, array_length, period)
% Sinusoidal - Computes a sinusoidal signal based on input parameters.
%   This function computes a sinusoidal signal based on the input array,
%   its length, and the period of the sinusoidal signal.
%
%   Inputs:
%     - input_array: Input array containing data points
%     - array_length: Length of the input array
%     - period: Period of the sinusoidal signal
%
%   Output:
%     - sinusoidal_signal: Array containing the computed sinusoidal signal

    % Calculate coefficients for sinusoidal signal
    coefficients = calculate_coefficients(input_array, array_length, period);
    b1 = coefficients(1);
    b2 = coefficients(2);
    
    % Initialize sinusoidal signal array
    sinusoidal_signal = zeros(1, array_length);
    
    % Compute sinusoidal signal for each time point
    for t = 1:array_length
        sinusoidal_signal(1, t) = b1 * cos(2 * pi * t / period) + b2 * sin(2 * pi * t / period);
    end
end

function coefficients = calculate_coefficients(input_array, array_length, period)
% calculate_coefficients - Calculates coefficients for the sinusoidal signal.
%   This function calculates the coefficients required to generate a sinusoidal
%   signal from the input array and its period.
%
%   Inputs:
%     - input_array: Input array containing data points
%     - array_length: Length of the input array
%     - period: Period of the sinusoidal signal
%
%   Output:
%     - coefficients: Array containing the calculated coefficients

    b1 = 0;
    b2 = 0;
    
    % Compute coefficients for the sinusoidal signal
    for t = 1:array_length
        b1 = b1 + input_array(t, 1) * cos(2 * pi * t / period);
        b2 = b2 + input_array(t, 1) * sin(2 * pi * t / period);
    end
    
    % Normalize coefficients
    coefficients(1) = b1 * 2 / array_length;
    coefficients(2) = b2 * 2 / array_length;
end
