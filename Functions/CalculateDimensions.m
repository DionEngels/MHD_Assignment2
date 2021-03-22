function [sym, res] = CalculateDimensions(sym, res, fixed, phys, coil, settings)
% CalculateDimensions Calculates the mayor radius, plasma volume, total cost, the total magnet thickness and the normal
% stress to which the structural material in the TF field coils is subjected
syms a

if settings.roger
    syms c
    res.epsilon = 0.4;
    sym.R_0 = 10 / sym.a;
    F_m = 2 * pi * sym.R_0 * coil.B_max^2 / phys.mu_0 * (2 * sym.a + 2 * res.b + c) / 2;
    F_t = coil.sigma_max * 4 * pi * sym.R_0 * c;
    sym.c = solve(F_m == F_t, c);
    
    cost_function = 2 * pi^2 * sym.R_0 * ((sym.a + res.b + sym.c)^2 - sym.a^2) / (fixed.P_e/1e6);
    
    res.(char(settings.var)) = fminbnd(matlabFunction(cost_function), settings.values(1), settings.values(end));
else
    sym.R_0 = (4/5 * phys.E_f * fixed.P_e)/(4 * pi^2 * fixed.eta_e * (phys.E_f + phys.E_Li)*sqrt((1 + res.kappa^2) / 2) * fixed.P_w * sym.a);
    
    xi = coil.B_max^2 / (4 * phys.mu_0 * phys.E_s * coil.e_crit);
    sigma = phys.E_s * coil.e_crit;
    nu = coil.B_max / (phys.mu_0 * coil.j_max) * (1 - coil.E/phys.E_s);
    
    sym.c = (nu + 2 * xi * (res.kappa * sym.a + res.b))/(1-xi);
    
    volumeperwatt = sqrt(2 / (1 + res.kappa^2)) * 4/5 * phys.E_f * ((res.kappa * sym.a + res.b + sym.c)*(sym.a + res.b + sym.c) - res.kappa * sym.a^2)/(2 * fixed.eta_e * (phys.E_f + phys.E_Li) * sym.a * fixed.P_w);
    
    try
        res.(char(settings.var)) = double(vpasolve(diff(volumeperwatt,a)==0,[settings.values(1) settings.values(end)]));
    end
end

sym.B_0 = coil.B_max * (sym.R_0 - res.a - res.b) / sym.R_0;


% Numerical minimization of the cost to find the cost optimizing value of a, while
% maximizing the usage of magnets


% fill in other parameters
res.c = double(subs(sym.c, a, res.a));
res.R_0 = double(subs(sym.R_0, a, res.a));
res.B_0 = double(subs(sym.B_0, a, res.a));
res.epsilon = res.a / res.R_0;

%% try to make double
names = fieldnames(res);
for k=1:numel(names)
    try
        res.(names{k}) = double(res.(names{k}));
    end
end

end