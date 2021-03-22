function [sym, res] = CalculateDimensions(sym, res, fixed, phys, coil, settings)
% CalculateDimensions Calculates the mayor radius, plasma volume, total cost, the total magnet thickness and the normal
% stress to which the structural material in the TF field coils is subjected
syms a c

if settings.roger
%     res.epsilon = 0.4;
%     res.R_0 = res.a / res.epsilon;
%     F_m = 2 * pi * res.R_0 * coil.B_max^2 / phys.mu_0 * (2 * res.a + 2 * res.b + c) / 2;
%     F_t = coil.sigma_max * 4 * pi * res.R_0 * c;
%     res.c = solve(F_m == F_t, c);
%     
%     cost_function = 2 * pi^2 * res.R_0 * ((res.a + res.b + res.c)^2 - res.a^2) / (fixed.P_e/1e6);
%     
%     figure;
%     fplot(cost_function, [1, 3], 'DisplayName', '13 T');
%     ylim([0, 4]);
%     
%     res.aval = fminbnd(matlabFunction(cost_function), 1, 4);
else
    sym.R_0 = (4/5 * phys.E_f * fixed.P_e)/(4 * pi^2 * fixed.eta_e * (phys.E_f + phys.E_Li)*sqrt((1 + res.kappa^2) / 2) * fixed.P_w * sym.a);

    xi = coil.B_max^2 / (4 * phys.mu_0 * phys.E_s * coil.e_crit);
    sigma = phys.E_s * coil.e_crit;
    nu = coil.B_max / (phys.mu_0 * coil.j_max) * (1 - coil.E/phys.E_s);
    
    sym.c = (nu + 2 * xi * (res.kappa * sym.a + res.b))/(1-xi);
    
    volumeperwatt = sqrt(2 / (1 + res.kappa^2)) * 4/5 * phys.E_f * ((res.kappa * sym.a + res.b + sym.c)*(sym.a + res.b + sym.c) - res.kappa * sym.a^2)/(2 * fixed.eta_e * (phys.E_f + phys.E_Li) * sym.a * fixed.P_w);

    res.a = double(vpasolve(diff(volumeperwatt,a)==0,[0 5]));
end

sym.B_0 = coil.B_max * (sym.R_0 - res.a - res.b) / sym.R_0;


% Numerical minimization of the cost to find the cost optimizing value of a, while
% maximizing the usage of magnets


% fill in other parameters
res.c = double(subs(sym.c, a, res.a));
res.R_0 = double(subs(sym.R_0, a, res.a));
res.B_0 = double(subs(sym.B_0, a, res.a));
end