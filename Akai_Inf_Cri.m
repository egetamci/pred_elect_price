
function [Aic_qp, ppaq] = Akai_Inf_Cri(p_rq, Y66q, N_s2q)
% PPA - Parameter Selection using Akaike Information Criterion (AIC) for α-stable PAR model
%   [Aic_qp, ppaq] = PPA(p_rq, Y66q, N_s2q)
%   This function performs parameter selection for a α-stable PAR model
%   using the Akaike Information Criterion (AIC).
%
%   Inputs:
%     - p_rq: Number of parameters
%     - Y66q: Input data
%     - N_s2q: Data length
%
%   Outputs:
%     - Aic_qp: Vector containing AIC values for different parameter selections
%     - ppaq: The value of ppq corresponding to the minimum AIC value
%
    % Initialize minimum AIC value
    min_pp = 10^(10);
    % Initialize variable to store the value of ppq corresponding to minimum AIC
    enk = 0;
    % Set the maximum value of ppq
    px = 15;
    
  % Loop through different values of ppq
    for ppq = 1:px
        % Initialize coefficient matrix
        aaq = -ones(ppq + 1, p_rq);
        % Calculate coefficients using Co_eff function
        coeffq = Co_eff(ppq, Y66q, N_s2q, p_rq);
        % Update the coefficient matrix
        aaq(2:end, :) = coeffq';
        % Calculate noise using APP_noise function
        Noise_new77q = APP_noise(ppq, aaq, Y66q, N_s2q, p_rq);
        % Estimate parameters using stblfit function
        alfa_pq = stblfit(Noise_new77q);
        % Calculate AIC using Aic_q function
        Aic_qp(ppq, 1) = Aic_q(Noise_new77q, alfa_pq(1, 1), 0, 1, 0, p_rq, ppq);
        
        % Update minimum AIC value and corresponding ppq value if necessary
        if (Aic_qp(ppq, 1) < min_pp)
            min_pp = Aic_qp(ppq, 1);
            enk = ppq;
        end
    end
    % Store the value of ppq corresponding to minimum AIC
    ppaq = enk;
end

function AIC = Aic_q(Noise_a, alfa_a, beta_a, gamma_a, Mu_a, peri_a, pp_a)
    % This function calculates the Akaike Information Criterion (AIC) value
    % Inputs:
    %   - Noise_a: Noisy data
    %   - alfa_a, beta_a, gamma_a, Mu_a: Parameters of the α-stable distribution
    %   - peri_a: Length of the time period
    %   - pp_a: Number of model parameters
    % Output:
    %   - AIC: Akaike Information Criterion value

    % Calculate the probability density function of the α-stable distribution
    alfa_stab = stblpdf(Noise_a, alfa_a, beta_a, gamma_a, Mu_a);
    % Take the natural logarithm of the probability density function
    alfa_stabg = log(alfa_stab);
    % Calculate the AIC value
    AIC = -2 * (sum(alfa_stabg)) + 2 * peri_a * (pp_a + 1);
end

function coeff5 = Co_eff(ppc, Yc, N_s2c, p_rc)
% Co_eff - Calculate Coefficients for AIC Calculation
%   Co_eff function calculates the coefficients matrix required for AIC calculation.
%
%   Inputs:
%     - ppc: Number of parameters
%     - Yc: Input data
%     - N_s2c: Data length
%     - p_rc: Number of iterations
%
%   Output:
%     - coeff5: Coefficients matrix
%
    % Initialize matrices for covariance and R values
    Lg_m1c = zeros(ppc, ppc);
    Rg_m1c = zeros(ppc, 1);
    
    % Loop through iterations
    for ii = 1:p_rc
        % Compute covariance matrix using cov_M function
        Lg_m1c(:,:) = cov_M(ii, ppc, Yc, N_s2c, p_rc);
        % Compute R vector using R_M function
        Rg_m1c(:, 1) = R_M(ppc, Yc, N_s2c, p_rc, ii);
        % Compute coefficients using matrix inversion
        coeff5(ii, :) = inv(Lg_m1c) * Rg_m1c;
    end
end


function cov = cov_M(ss1, pp1, X33, N33, T33)
% cov_M - Calculate Covariance Matrix
%   cov_M function computes the covariance matrix.
%
%   Inputs:
%     - ss1: Some parameter
%     - pp1: Number of parameters
%     - X33, N33, T33: Other parameters
%
%   Output:
%     - cov: Covariance matrix
%
    % Initialize covariance matrix
    cov = zeros(pp1, pp1);
    
    % Loop through each element of the covariance matrix
    for i3 = 1:pp1
        for j3 = 1:pp1
            % Compute each element using Lamda_M function
            cov(i3, j3) = Lamda_M(X33, T33, ss1 - j3, i3 - j3, N33);
        end
    end
end



