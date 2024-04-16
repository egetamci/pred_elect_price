function L_M1 = Covariance(X1, p)
% Computes the covariance sequence of a given signal.
% Inputs:
%   - X1: Input signal vector
%   - p: Number of lags for computing covariance
% Output:
%   - L_M1: Covariance sequence
x_l = length(X1);
Top_X = sum(abs(X1));
for i2 = 1:p
L_T1 = 0;
for k3 = i2 + 1:x_l
L_T1_aa = X1(k3) * sign(X1(k3 - i2));
L_T1 = L_T1 + L_T1_aa;
end
L_M1(i2) = L_T1;
end
L_M1 = L_M1 ./ Top_X;
end
