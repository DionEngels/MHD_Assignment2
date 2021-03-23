close all;
clear;
clc;
addpath('Functions');

%%%% IMPORTANT: convert every parameter and formula to SI units before
%%%% using it in this program %%%%%%

%% Settings

settings.use_rebco = 0;               % Use REBCO magnets or not
settings.roger = 0;                   % 1 = use sigma constraint on magnets, 0 = use epsilon_irr - j_max constarint on magnets
settings.variable = 1;                % 1 = a, 2 = P_e, 3 = P_w, 4 = H
settings.plotvolumeperwatt = 0;       % 1 = yes, 0 = no

%% Physical constants

phys.mu_0 = 1.256e-6;                 % H/m
phys.e = 1.60217649e-19;              % C
phys.k = 1.3806504e-23;               % J*C/K
phys.E_f = 17.6e6 * phys.e;           % J energy per fusion reaction
phys.E_Li = 4.8e6 * phys.e;           % J energy per neutron produced via Li-6 breeding in the blanket
phys.E_s = 210e9;                     % Pa Youngs modulus of the structural material (steel)
phys.sigma = 300e6;                   % Pa Tensile strength of the structural material
phys.sigma_v = 3e-22;                 % m3/s (at 15 keV)
phys.A = 2.5;                         % - average atomic mass for D-T

%% Materials data from the superconducting magnets

coil.Nb.B_max = 13;                   % T
coil.Re.B_max = 19;                   % T
coil.Nb.j_max = 1000e6;               % A/m^2
coil.Re.j_max = 2500e6;               % A/m^2
coil.Nb.e_crit =  0.2/100;            % - (not % !)
coil.Re.e_crit =  0.62/100;           % - (not % !)
coil.Nb.E =  150e9;                   % Pa
coil.Re.E =  172e9;                   % Pa

coil.Nb.sigma_max = 300e6;            % m3/s (at 15 keV)
coil.Re.sigma_max = 300e6;            % m3/s (at 15 keV)

coil.High.B_max = 20;                 % T For plotting
coil.Low.B_max = 6;                   % T For plotting

%% Predetermined parameters/results that aren't optimized

res.kappa = 1.7;                    % - (Stabilization limit of the n=0 vertical displacement mode)
res.b = 1.2;                        % m (blanket thickness based on nuclear engineering constraints)
res.T = 15e3 * phys.e / phys.k;     % K (Temperature at which the fusion power is optimized for a fixed pressure)
res.eta_e = 0.4;                    % - thermal conversion efficiency

%% Parameter to vary

switch settings.variable
    
    case 1
        
        syms a                 
        res.a = a;
        settings.costoptimizea = 0;
        
        res.P_e = 1000e6;
        res.P_w = 4e6;
        res.H = 1;
        
        settings.var = a;
        settings.values = linspace(1,1.7,100);  % Range of a values to plot and the resolution
        
    case 2
        
        syms P_e a
        
        res.a = a;
        settings.costoptimizea = 1;
        
        res.P_e = P_e;
        
        res.P_w = 4e6;
        res.H = 1;
        
        settings.var = P_e;
        settings.values = linspace(800,2000,100)*1e6;  % Range of a values to plot and the resolution
        
    case 3
        
        syms P_w a
        
        res.a = a;
        settings.costoptimizea = 1;
        
        res.P_w = P_w;
        res.P_e = 1000e6;
        res.H = 1;
        
        settings.var = P_w;
        settings.values = linspace(2,5,100)*1e6;  % Range of a values to plot and the resolution
        
    case 4
        
        syms H a
        
        res.a = a;
        settings.costoptimizea = 1;
        
        res.H = H;
        
        res.P_e = 1000e6;
        res.P_w = 4e6;
        
        settings.vartocostoptimize = a;
        
        settings.var = H;
        settings.values = linspace(0.5,5,100);  % Range of a values to plot and the resolution
        
    otherwise
        
        disp("You selected invalid variable case");
end

%% Volume and cost of the device from the wall loading and power output contraints (engineering)
if ~settings.use_rebco
    
    [res] = CalculateDimensions(res, phys, coil.Nb, settings);
    
else
    
    [res] = CalculateDimensions(res, phys, coil.Re, settings);
    
end

%% Resulting plasma parameters based on scaling laws, steady state ignition and the required power output
[res] = CalculatePlasma(res, phys, settings);

%% Checking if the required plasma parameters are feasible through the normalized limits (plasma physics)
[lim] = ConstraintsCheck(res, phys);

%% Plot
ConstraintsPlot(lim, res, settings);

%% Print Results
PrintResults(res);
