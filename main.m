close all; addpath('Functions');
%% Settings
settings.a_bounds = [0, 4];           % bounds for minor radius in plots
settings.constraint_plot = 1;         % 1 for show plots, 0 for calculate check constraints
settings.use_rebco = 0;               % Use REBCO magnets or not

%% Physical Constants
phys.mu_0 = 1.256e-6;                 % H/m
phys.e = 1.60217649e-19;              % C
phys.k = 1.3806504e-23;               % J*C/K
phys.E_f = 17.6e6 * phys.e;           % J energy per fusion reaction
phys.E_Li = 4.8e6 * phys.e;           % J energy per neutron produced via Li-6 breeding in the blanket
phys.sigma_sd = 1e-28;                % m2 Slowdown crosssection
phys.sigma_br = 950e-28;              % m2 Breeding crosssection
phys.sigma_v = 3e-22;                 % m3/s (at 15 keV)
phys.lambda_mfp = 0.02;               % m mean free path
phys.A = 2.5;                         % - average atomic mass for D-T
phys.H = 1;                           % - H mode enhancement factor (tau_E is multiplied with this factor)

%% Fixed parameters
fixed.P_e = 1000e6;                   % W power provided to the grid
fixed.P_w = 4e6;                      % W/m2 max neutron wall loading
fixed.eta_e = 0.4;                    % - thermal conversion efficiency
%fixed.E_s = 210e9;                    % Pa Youngs modulus of the structural material (steel)

%% Materials data from the superconducting magnets

coil.Nb.B_max = 13;                   % T
coil.Re.B_max = 19;                   % T
coil.Nb.j_max = 1000e6;               % A/m^2
coil.Re.j_max = 2000e6;               % A/m^2
coil.Nb.e_crit =  0.8/100;            % - (not % !)
coil.Re.e_crit =  0.55/100;           % - (not % !)
coil.Nb.E =  120e9;                   % Pa
coil.Re.E =  170e9;                   % Pa
coil.Nb.sigma_max = 300e6;            % m3/s (at 15 keV)
coil.Re.sigma_max = 300e6;            % m3/s (at 15 keV)

coil.High.B_max = 20;                 % T For plotting
coil.Low.B_max = 6;                   % T For plotting

%% Predetermined parameters that aren't optimized
res.kappa = 1.7;                    % - (Stabilization limit of the n=0 vertical displacement mode)
res.b = 1.2;                        % m (blanket thickness based on nuclear engineering constraints)
res.T = 15e3 * phys.e / phys.k;     % K (Temperature at which the fusion power is optimized for a fixed pressure)
res.epsilon = 0.4;                  % inverse aspect ratio
res.P_fusion = fixed.P_e / fixed.eta_e; % Fusion power

%% Desired parameter
syms a;                               % m
res.R_0 = a / res.epsilon;            % m
%% Required parameter
syms a;                 % m

%% Unknowns
R_0 = 10 / a;           % m

%% Find c
if ~settings.use_rebco
    c = CoilForceBalance(a, R_0, phys, coil.Nb, res);
else
    c = CoilForceBalance(a, R_0, phys, coil.Re, res);
end

%% Plot and find a using Cost function
figure;
title("Cost function again minor radius")
hold on
CostFunction(a, res.b, c, R_0, fixed, coil.Low, settings, '6 T');
if ~settings.use_rebco
    res.a = CostFunction(a, res.b, c, R_0, fixed, coil.Nb, settings, '13 T');
else
    res.a = CostFunction(a, res.b, c, R_0, fixed, coil.Re, settings, '19 T');
end
CostFunction(a, res.b, c, R_0, fixed, coil.High, settings, '20 T');
legend
hold off

%% Evaluate design
if ~settings.use_rebco
    res = EvaluateDesign(res, coil.Nb, c, R_0, settings);
else
    res = EvaluateDesign(res, coil.Re, c, R_0, settings);
end

clear c ans R_0

%% Plasma requirements
res = CalculatePlasmaParameters(fixed, phys, res);

%% Constraints
ConstraintsCheck(res, settings)

%% Print found parameters
PrintResults(res)
