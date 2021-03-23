function [res] = CalculateDimensions(res, phys, coil, settings)

% CalculateDimensions Calculates the mayor radius, plasma volume, total cost, the total magnet thickness and the normal
% stress to which the structural material in the TF field coils is subjected

if settings.roger
    
    res.R_0 = (4/5 * phys.E_f * res.P_e)/(4 * pi^2 * res.eta_e * (phys.E_f + phys.E_Li)*sqrt((1 + res.kappa^2) / 2) * res.P_w * res.a);
    
    xi = coil.B_max^2 / (4 * phys.mu_0 * phys.sigma);
    
    res.c = 2 * xi * (res.kappa * res.a + res.b)/(1-xi); 
    
    res.B_0 = coil.B_max * (res.R_0 - res.a - res.b) / res.R_0;
    
    res.volumeperwatt = sqrt(2 / (1 + res.kappa^2)) * 4/5 * phys.E_f * ((res.kappa * res.a + res.b + res.c)*(res.a + res.b + res.c) - res.kappa * res.a^2)/(2 * res.eta_e * (phys.E_f + phys.E_Li) * res.a * res.P_w);
    
    if settings.costoptimizea
        syms a
        res.a = solve(diff(res.volumeperwatt,res.a)==0,a);
        res.a = res.a(double(res.a)>=0);
        res.c = subs(res.c, a, res.a);
        res.R_0 = subs(res.R_0, a, res.a);
        res.B_0 = subs(res.B_0, a, res.a);
        res.volumeperwatt = subs(res.volumeperwatt, a, res.a);
        res.epsilon = res.a / res.R_0;
    end
    
else
    
    res.R_0 = (4/5 * phys.E_f * res.P_e)/(4 * pi^2 * res.eta_e * (phys.E_f + phys.E_Li)*sqrt((1 + res.kappa^2) / 2) * res.P_w * res.a);
    
    xi = coil.B_max^2 / (4 * phys.mu_0 * phys.E_s * coil.e_crit);
    
    nu = coil.B_max / (phys.mu_0 * coil.j_max) * (1 - coil.E/phys.E_s);
    
    res.c = (nu + 2 * xi * (res.kappa * res.a + res.b))/(1-xi); 
    
    res.B_0 = coil.B_max * (res.R_0 - res.a - res.b) / res.R_0;
    
    res.volumeperwatt = sqrt(2 / (1 + res.kappa^2)) * 4/5 * phys.E_f * ((res.kappa * res.a + res.b + res.c)*(res.a + res.b + res.c) - res.kappa * res.a^2)/(2 * res.eta_e * (phys.E_f + phys.E_Li) * res.a * res.P_w);
    
    if settings.costoptimizea
        syms a
        res.a = solve(diff(res.volumeperwatt,res.a)==0,a);
        res.a = res.a(double(res.a)>=0);
        res.c = subs(res.c, a, res.a);
        res.R_0 = subs(res.R_0, a, res.a);
        res.B_0 = subs(res.B_0, a, res.a);
        res.volumeperwatt = subs(res.volumeperwatt, a, res.a);
        res.epsilon = res.a / res.R_0;
    end
end



%% try to make numerical values (doubles)

names = fieldnames(res);
for k=1:numel(names)
    try
        res.(names{k}) = double(res.(names{k}));
    end
end

end