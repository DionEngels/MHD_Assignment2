function [res] = CalculatePlasmaParameters(fixed, phys, res)
%CALCULATECURRENT Calculates plasma parameters such as I_p, p, n, tau_e
%% p
res.p = sqrt(16*fixed.P_e*(phys.k*res.T)^2/(fixed.eta_e*(phys.E_f+phys.E_Li)*res.V_p*phys.sigma_v));
%% beta
res.beta = res.p / (res.B_0^2 / (2 * phys.mu_0));
%% n
res.n = res.p / (2 * res.T * phys.k);

%% I_p
syms I_p
I_p = solve(1.2 == 0.082 * I_p * (res.P_fusion / 1e6)^(-0.5) * res.R_0^(1.6) * res.kappa^(-0.2) * res.B_0^(0.15) * phys.A^(0.5), I_p);

try
    I_p = double(I_p);
catch
    I_p = I_p;
end
res.I_p = I_p * 1e6; % from MA to A

%% tau_e
res.tau_e = 3/2 * res.p / (1/4* res.n^2 * phys.sigma_v * 1/5 * phys.E_f);
end

