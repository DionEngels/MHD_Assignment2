close all
%% Settings
settings_a_bounds = [0, 4];     % bounds for minor radius in plots
settings_constraint_plot = 0;   % 1 for show plots, 0 for calculate check constraints

%% Given parameters
P_e = 1000;             % MW
P_w = 4;                % MW/m2
syms B_max;             % T
B_max_val = 13;         % T
sigma_max = 300e6;      % Pa
reactivity_dt = 3e-22;  % m3/s
sigma_sd = 1e-28;       % m2
sigma_br = 950e-28;     % m2
eta_e = 0.4;            % -
A = 1;                  % kg
lambda_mfp = 0.02;       % m

%% Required parameters
syms a;                 % m
R_0 = 10 / a;           % m
epsilon = 0.1;          % -
b = 1.2;                % m
syms c;                 % m
syms V_p;               % m3
syms A_p;               % m2
syms B_0;               % T
syms I_p;               % MA
syms p;                 % Pa
syms T;                 % eV
syms n;                 % /m3
syms tau_e;             % s
syms beta;              % -
kappa = 1.7;             % -
syms P_fusion;          % MW
P_fusion = P_e / eta_e;

%% Physical values
mu_0 = 1.256637e-6;     % H/m

%% Find c
c = CoilForceBalance(R_0, B_max, mu_0, a, b, sigma_max);

%% Cost function
cost_function = 2 * pi^2 * R_0 * ((a + b + c)^2 - a^2) / P_e;

%% plot cost versus a
figure;
fplot(subs(cost_function, B_max, 6), settings_a_bounds);
hold on
fplot(subs(cost_function, B_max, 13), settings_a_bounds);
fplot(subs(cost_function, B_max, 20), settings_a_bounds);
ylim([0, 4]);

%% Evaluate design
cost_function = subs(cost_function, B_max, B_max_val);
c = subs(c, B_max, B_max_val);
a_val = 2;
c = double(subs(c, a, a_val));
R_0 = double(subs(R_0, a, a_val));
B_max = B_max_val;
if ~settings_constraint_plot
    a = a_val;
end
B_0 = ((a + b + c/2) / 2) * B_max / R_0;

%% Plasma requirements
p = 7e5;
T = 15e3;
n = 1.4e20;
beta = 0.08;
tau_e = 1.2;

%% Confinement
if ~settings_constraint_plot
    I_p = double(solve(1.2 == 0.082 * I_p * P_fusion^(-0.5) * R_0^(1.6) * kappa^(-0.2) * B_0^(0.15) * A^(0.5), I_p)); % MA
else
    I_p = solve(1.2 == 0.082 * I_p * P_fusion^(-0.5) * R_0^(1.6) * kappa^(-0.2) * B_0^(0.15) * A^(0.5), I_p); % MA
end

%% Constraints
q = (5 * a^2 * kappa * B_0) / (R_0 * I_p);
f_b = 1.3 * kappa^(1/4) * beta * (5 * a^2 * kappa * B_0) / (R_0 * I_p) / epsilon^(1/2);

limit_kappa = 2;
limit_q = 2;
limit_troyon = 0.028 * I_p / a / B_0;
limit_greenwald = I_p / (pi * a^2);
limit_bootstrap = 0.8;
%% Plot constraints or calculate constraints
if settings_constraint_plot
    figure;
    hold on
    fplot(f_b / limit_bootstrap, settings_a_bounds);
    fplot(n / limit_greenwald, settings_a_bounds);
    fplot(beta / limit_troyon, settings_a_bounds);
    fplot(q / limit_q, settings_a_bounds);
    ylim([0, 2])
else
    LogicalStr = {'false', 'true'};
    
    constraint_kappa = kappa < limit_kappa;
    constraint_safety = q > limit_q;
    constraint_troyon = beta < limit_troyon;
    constraint_greenwald = n < limit_greenwald;
    constraint_bootstrap = f_b > limit_bootstrap;
    
    n_constraints_passed = constraint_kappa + constraint_safety + constraint_troyon + constraint_greenwald + constraint_bootstrap;
    fprintf('Kappa constraint is fulfilled: %s\n', LogicalStr{constraint_kappa + 1})
    fprintf('Safety factor constraint is fulfilled: %s\n', LogicalStr{constraint_safety + 1})
    fprintf('Troyon limit constraint is fulfilled: %s\n', LogicalStr{constraint_troyon + 1})
    fprintf('Greenwald limit constraint is fulfilled: %s\n', LogicalStr{constraint_greenwald + 1})
    fprintf('Bootstrap current constraint is fulfilled: %s\n', LogicalStr{constraint_bootstrap + 1})
    fprintf('Result: %u passed, %u failed\n', n_constraints_passed, 5 - n_constraints_passed)
end
