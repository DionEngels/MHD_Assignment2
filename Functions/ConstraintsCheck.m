function [] = ConstraintsCheck(sym, res, phys, settings)
% ConstraintsCheck Calculates the normalized ratio of the plasma parameter with prominent
% limits. There is compliance if all limits are < 1
syms a

% External kink safety factor
q_star_sym = 2 * pi * sym.a^2 * sym.B_0 / (phys.mu_0 * sym.R_0 * sym.I_p) * ((1 + res.kappa^2) / 2);
q_star_val = 2 * pi * res.a^2 * res.B_0 / (phys.mu_0 * res.R_0 * res.I_p) * ((1 + res.kappa^2) / 2);

lim.sym.q = 2/q_star_sym;
lim.res.q = 2/q_star_val;

% Troyon beta limit

beta_val = res.p / (res.B_0^2 / (2 * phys.mu_0));
beta_sym = sym.p / (sym.B_0^2 / (2 * phys.mu_0));

beta_troyon_sym = 2.8 / 100 * (sym.I_p * 1e-6 / (sym.a * sym.B_0));
beta_troyon_val = 2.8 / 100 * (res.I_p * 1e-6 / (res.a * res.B_0));

lim.sym.troyon = beta_sym/beta_troyon_sym;
lim.res.troyon = beta_val/beta_troyon_val;

% Greenwald density limit

lim.res.greenwald = res.n / (res.I_p * 1e-6 / (pi * res.a^2) * 1e20);
lim.sym.greenwald = sym.n / (sym.I_p * 1e-6 / (pi * sym.a^2) * 1e20);

% Required neoclassical bootstrap current fraction (profile assumptions for
% every plasma parameter, see pages 525-526 of FB)

lim.res.bootstrap = 0.75/ (4/3 * res.kappa^(1.4) * 2.8/100 * q_star_val / (res.a / res.R_0)^0.5);
lim.sym.bootstrap = 0.75/ (4/3 * res.kappa^(1.4) * 2.8/100 * q_star_sym / (sym.a / sym.R_0)^0.5);

%% Plot and print
ConstraintsPlot(lim.sym, settings);
ConstraintsPrint(lim.res);

end