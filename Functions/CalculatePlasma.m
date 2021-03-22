function [res] = CalculatePlasma(res, fixed, phys)
% CalculatePlasma Calculates p, n, tau_E and I_p of the plasma based on the required
% power output, steady state ignition criterium and ELMy H mode scaling law and the predetermined volume of the plasma

res.V_p =  2 * pi^2 * res.R_0 * res.kappa * res.a^2;

res.p = sqrt(16 * fixed.P_e * (phys.k * res.T)^2 / (fixed.eta_e * (phys.E_f + phys.E_Li) * res.V_p * phys.sigma_v));
res.n = res.p/(2 * phys.k * res.T);

res.tau_e = 3/2 * res.p / (1/4 * res.n^2 * phys.sigma_v * 1/5 * phys.E_f);

P_alpha = 1/4 * res.n^2 * phys.sigma_v * 1/5 * phys.E_f * res.V_p;
res.I_p = ((res.tau_e * (P_alpha * 1e-6)^(0.69)) / (0.145 * phys.H * res.R_0^(1.39) * res.a^(0.58) * res.kappa^(0.78) * (res.n*1e-20)^(0.41) * res.B_0^(0.15) * phys.A^(0.19)))^(1/0.93)*1e6;

end