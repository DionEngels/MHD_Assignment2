function [] = ConstraintsCheck(res, phys, settings)
% ConstraintsCheck Calculates the normalized ratio of the plasma parameter with prominent
% limits. There is compliance if all limits are < 1
syms a

% External kink safety factor
q_star = 2 * pi * res.a^2 * res.B_0 / (phys.mu_0 * res.R_0 * res.I_p) * ((1 + res.kappa^2) / 2);

lim.q = 2/q_star;

% Troyon beta limit

beta = res.p / (res.B_0^2 / (2 * phys.mu_0));

beta_troyon = 2.8 / 100 * (res.I_p * 1e-6 / (res.a * res.B_0));

lim.troyon = beta/beta_troyon;

% Greenwald density limit

n_greenwald = res.I_p * 1e-6 / (pi * res.a^2) * 1e20;

lim.greenwald = res.n / n_greenwald;

% Required neoclassical bootstrap current fraction (profile assumptions for
% every plasma parameter, see pages 525-526 of FB)

f_b = 4/3 * res.kappa^(1.4) * 2.8/100 * q_star / (res.a / res.R_0)^0.5;

lim.bootstrap = 0.75/f_b;

%% Plot and print
ConstraintsPlot(lim, settings);
ConstraintsPrint(res, lim);

end