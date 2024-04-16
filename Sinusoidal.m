function num_fonk = Sinusoidal(array2, l_n2, p_r2)
% This function computes a sinusoidal signal based on input parameters.
% Inputs:
%   - array2: Input array
%   - l_n2: Length of the array
%   - p_r2: Period of the sinusoidal signal
% Output:
%   - num_fonk1: Sinusoidal signal array
coeffM = calculate_coeff_sinu(array2, l_n2, p_r2);
b11 = coeffM(1);
b22 = coeffM(2);
for t2 = 1:l_n2
num_fonk(1, t2) = b11 * cos(2 * pi * t2 / p_r2) + b22 * sin(2 * pi * t2 / p_r2);
end
end
function Coef = calculate_coeff_sinu(array1, l_n1, p_r1)
b11 = 0;
b22 = 0;
for tt = 1:l_n1
b11 = b11 + array1(tt, 1) * cos(2 * pi * tt / p_r1);
b22 = b22 + array1(tt, 1) * sin(2 * pi * tt / p_r1);
end
Coef(1) = b11 * 2 / l_n1;
Coef(2) = b22 * 2 / l_n1;
end
