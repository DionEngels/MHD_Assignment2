function [I_p] = CalculateCurrent(fixed, res)
%CALCULATECURRENT Calculates the I_p
%   Detailed explanation goes here
syms I_p
P_fusion = fixed.P_e / (1e6) / fixed.eta_e;
I_p = solve(1.2 == 0.082 * I_p * P_fusion^(-0.5) * res.R_0^(1.6) * res.kappa^(-0.2) * res.B_0^(0.15) * fixed.A^(0.5), I_p);

try
    I_p = double(I_p);
catch
    I_p = I_p;
end

