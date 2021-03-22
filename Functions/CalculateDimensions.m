function [res] = CalculateDimensions(res, fixed, phys, coil)
% CalculateDimensions Calculates the mayor radius, plasma volume, total cost, the total magnet thickness and the normal
% stress to which the structural material in the TF field coils is subjected
syms a

res.R_0 = (4/5 * phys.E_f * fixed.P_e)/(4 * pi^2 * fixed.eta_e * (phys.E_f + phys.E_Li)*sqrt((1 + res.kappa^2) / 2) * fixed.P_w * res.a);

res.B_0 = coil.B_max * (res.R_0 - res.a - res.b) / res.R_0;

xi = coil.B_max^2 / (4 * phys.mu_0 * phys.E_s * coil.e_crit);
sigma = phys.E_s * coil.e_crit;
nu = coil.B_max / (phys.mu_0 * coil.j_max) * (1 - coil.E/phys.E_s);

res.c = (nu + 2 * xi * (res.kappa * res.a + res.b))/(1-xi);

volumeperwatt = sqrt(2 / (1 + res.kappa^2)) * 4/5 * phys.E_f * ((res.kappa * res.a + res.b + res.c)*(res.a + res.b + res.c) - res.kappa * res.a^2)/(2 * fixed.eta_e * (phys.E_f + phys.E_Li) * res.a * fixed.P_w);

% Numerical minimization of the cost to find the cost optimizing value of a, while
% maximizing the usage of magnets

res.aval = double(vpasolve(diff(volumeperwatt,a)==0,[0 5]));

% fill in other parameters
res.cval = double(subs(res.c, a, res.aval));
res.R_0val = double(subs(res.R_0, a, res.aval));
res.B_0val = double(subs(res.B_0, a, res.aval));
end