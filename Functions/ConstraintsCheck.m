function [lim] = ConstraintsCheck(res, phys)
% ConstraintsCheck Calculates the normalized ratio of the plasma parameter with prominent
% limits. There is compliance if all limits are < 1
syms a

% External kink safety factor
q_star_val = 2 * pi * res.a^2 * res.B_0 / (phys.mu_0 * res.R_0 * res.I_p) * ((1 + res.kappa^2) / 2);

lim.q = 2/q_star_val;

% Troyon beta limit

beta_val = res.p / (res.B_0^2 / (2 * phys.mu_0));

beta_troyon_val = 2.8 / 100 * (res.I_p * 1e-6 / (res.a * res.B_0));

lim.troyon = beta_val/beta_troyon_val;

% Greenwald density limit

lim.greenwald = res.n / (res.I_p * 1e-6 / (pi * res.a^2) * 1e20);

% Required neoclassical bootstrap current fraction (profile assumptions for
% every plasma parameter, see pages 525-526 of FB)

lim.bootstrap = 0.75/ (4/3 * res.kappa^(1.4) * 2.8/100 * q_star_val / (res.a / res.R_0)^0.5);

end