close all
%% Settings
settings.a_bounds = [0, 4];           % bounds for minor radius in plots
settings.constraint_plot = 1;         % 1 for show plots, 0 for calculate check constraints
settings.use_rebco = 0;               % Use REBCO magnets or not

%% Physical Constants
phys.mu_0 = 1.256e-6;                 % H/m
phys.e = 1.60217649e-19;              % C
phys.k = 1.3806504e-23;               % J*C/K

%% Fixed parameters
fixed.P_e = 1000e6;                   % W power provided to the grid
fixed.P_w = 4e6;                      % W/m2 max neutron wall loading
fixed.sigma_max = 300e6;              % m3/s (at 15 keV)
fixed.sigma_sd = 1e-28;               % m2 Slowdown crosssection
fixed.sigma_br = 950e-28;             % m2 Breeding crosssection
fixed.eta_e = 0.4;                    % - thermal conversion efficiency
fixed.A = 2.5;                        % - average atomic mass for D-T
fixed.H = 1;                          % - H mode enhancement factor (tau_E is multiplied with this factor)
fixed.E_f = 17.6e6 * phys.e;            % J energy per fusion reaction
fixed.E_Li = 4.8e6 * phys.e;            % J energy per neutron produced via Li-6 breeding in the blanket
fixed.E_s = 210e9;                    % Pa Youngs modulus of the structural material (steel)
fixed.lambda_mfp = 0.02;              % m mean free path

%% Materials data from the superconducting magnets

coil.Nb.B_max = 13;                   % T
coil.Re.B_max = 19;                   % T
coil.Nb.j_max = 1000e6;               % A/m^2
coil.Re.j_max = 2000e6;               % A/m^2
coil.Nb.e_crit =  0.8/100;            % - (not % !)
coil.Re.e_crit =  0.55/100;           % - (not % !)
coil.Nb.E =  120e9;                   % Pa
coil.Re.E =  170e9;                   % Pa

coil.High.B_max = 20;                 % T For plotting
coil.Low.B_max = 6;                   % T For plotting

%% Predetermined parameters that aren't optimized
res.kappa = 1.7;                    % - (Stabilization limit of the n=0 vertical displacement mode)
res.b = 1.2;                        % m (blanket thickness based on nuclear engineering constraints)
res.T = 15e3 * phys.e / phys.k;     % K (Temperature at which the fusion power is optimized for a fixed pressure)
res.epsilon = 0.4;                  % inverse aspect ratio

%% Desired parameter
syms a;                               % m

res.R_0 = a / res.epsilon;          % m
%% Required parameter
syms a;                 % m

%% Unknowns
R_0 = 10 / a;           % m

%% Find c
c = CoilForceBalance(a, R_0, phys, fixed, res);

%% Plot and find a using Cost function
figure;
hold on
CostFunction(a, res.b, c, R_0, fixed, coil.Low, settings);
if ~settings.use_rebco
    res.a = CostFunction(a, res.b, c, R_0, fixed, coil.Nb, settings);
else
    res.a = CostFunction(a, res.b, c, R_0, fixed, coil.Re, settings);
end
CostFunction(a, res.b, c, R_0, fixed, coil.High, settings);
hold off

%% Evaluate design
if settings.use_rebco
    res = EvaluateDesign(res, coil.Re, c, R_0, settings);
else
    res = EvaluateDesign(res, coil.Nb, c, R_0, settings);
end

%% Plasma requirements
p = 7e5;
n = 1.4e20;
beta = 0.08;
tau_e = 1.2;

%% Confinement

if ~settings.constraint_plot
    res.I_p = CalculateCurrent(fixed, res); % MA
else
    I_p = CalculateCurrent(fixed, res); % MA
end

%% Constraints
q = (5 * a^2 * kappa * B_0) / (R_0 * I_p);
f_b = 1.3 * kappa^(1/4) * beta * (5 * a^2 * kappa * B_0) / (R_0 * I_p) / epsilon^(1/2);

limit.kappa = 2;
limit.q = 2;
limit.troyon = 0.028 * I_p / a / B_0;
limit.greenwald = I_p / (pi * a^2);
limit.bootstrap = 0.8;
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
