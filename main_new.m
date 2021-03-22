close all

%%%% IMPORTANT: convert every parameter and formula to SI units before
%%%% using it in this program %%%%%%

%% Settings

as = linspace(1,1.7,100); % Range of a values to plot and the resolution
settings.use_rebco = 0;               % Use REBCO magnets or not

%% Physical constants

phys.mu_0 = 1.256e-6;                 % H/m
phys.e = 1.60217649e-19;              % C
phys.k = 1.3806504e-23;               % J*C/K
phys.E_f = 17.6e6 * phys.e;           % J energy per fusion reaction
phys.E_Li = 4.8e6 * phys.e;           % J energy per neutron produced via Li-6 breeding in the blanket
phys.E_s = 210e9;         % Pa Youngs modulus of the structural material (steel) 
phys.sigma_v = 3e-22;   % m3/s (at 15 keV)
phys.A = 2.5;                % - average atomic mass for D-T 
phys.H = 1;                  % - H mode enhancement factor (tau_E is multiplied with this factor)

%% Fixed parameters

fixed.P_e = 1000e6;        % W power provided to the grid
fixed.P_w = 4e6;           % W/m2 max neutron wall loading
fixed.eta_e = 0.4;            % - thermal conversion efficiency


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

%% Predetermined parameters that aren't optimized
res.kappa = 1.7;                    % - (Stabilization limit of the n=0 vertical displacement mode)
res.b = 1.2;                        % m (blanket thickness based on nuclear engineering constraints)
res.T = 15e3 * phys.e / phys.k;     % K (Temperature at which the fusion power is optimized for a fixed pressure)
res.P_fusion = fixed.P_e / fixed.eta_e; % Fusion power

%% Parameter to vary
syms a;                 % m
res.a = a;

%% Volume and cost of the device from the wall loading and power output contraints (engineering)
if ~settings.use_rebco
    res = dimensions(res, fixed, phys, coil.Nb);
else
    res = dimensions(res, fixed, phys, coil.Re);
end

%% Resulting plasma parameters based on scaling laws, steady state ignition and the required power output
res = plasmaparameters(res, fixed, phys);

%% Checking if the required plasma parameters are feasible through the normalized limits (plasma physics)

plasmalimits(as, res, phys);

%% Plot constraints or calculate constraints

   
%     figure(1)
%     hold on
%     plot(as,double(subs(bootstrap_lim_Nb,a,as)),'-r','DisplayName','0.75/f_b Nb3Sn');
%     plot(as,subs(greenwald_lim_Nb,a,as),'--r','DisplayName','n/n_G Nb3Sn');
%     plot(as,subs(troyon_lim_Nb,a,as),'-.r','DisplayName','{\beta}/{\beta_T} Nb3Sn');
%     plot(as,subs(q_lim_Nb,a,as),':r','DisplayName','2/q_{*} Nb3Sn');  
%     plot(as,subs(bootstrap_lim_Re,a,as),'-b','DisplayName','0.75/f_b ReBCuO');
%     plot(as,subs(greenwald_lim_Re,a,as),'--b','DisplayName','n/n_G ReBCuO');
%     plot(as,subs(troyon_lim_Re,a,as),'-.b','DisplayName','{\beta}/{\beta_T} ReBCuO');
%     plot(as,subs(q_lim_Re,a,as),':b','DisplayName','2/q_{*} ReBCuO');
%     legend
%     hold off
%     
%     figure(2)
%     as = linspace(1,3,100); % Range of a values to plot and the resolution
%     hold on
%     plot(as,subs(volperwatt_Nb/(10^(-6)),a,as),'-r','DisplayName','Volume per Watt Nb3Sn [m^3/(MW)]');
%     plot(as,subs(volperwatt_Re/(10^(-6)),a,as),'-b','DisplayName','Volume per Watt ReBCuO [m^3/(MW)]');
%     legend
%     hold off
    
    % To be reexamined !
    
    % LogicalStr = {'false', 'true'};
   
%     constraint_safety_Nb = q_lim_Nb < 1;
%     constraint_safety_Re = q_lim_Re < 1;
%     constraint_troyon_Nb = troyon_lim_Nb < 1;
%     constraint_troyon_Re = troyon_lim_Re < 1;
%     constraint_greenwald_Nb = greenwald_lim_Nb < 1;
%     constraint_greenwald_Nb = greenwald_lim_Re < 1;
%     constraint_bootstrap = f_b > limit_bootstrap;
%     
%     n_constraints_passed = constraint_kappa + constraint_safety + constraint_troyon + constraint_greenwald + constraint_bootstrap;
%     fprintf('Kappa constraint is fulfilled: %s\n', LogicalStr{constraint_kappa + 1})
%     fprintf('Safety factor constraint is fulfilled: %s\n', LogicalStr{constraint_safety + 1})
%     fprintf('Troyon limit constraint is fulfilled: %s\n', LogicalStr{constraint_troyon + 1})
%     fprintf('Greenwald limit constraint is fulfilled: %s\n', LogicalStr{constraint_greenwald + 1})
%     fprintf('Bootstrap current constraint is fulfilled: %s\n', LogicalStr{constraint_bootstrap + 1})
%     fprintf('Result: %u passed, %u failed\n', n_constraints_passed, 5 - n_constraints_passed)

function [res] = dimensions(res, fixed, phys, coil)
% COSTANDVOLUME Calculates the mayor radius, plasma volume, total cost, the total magnet thickness and the normal
% stress to which the structural material in the TF field coils is subjected
syms a

res.R_0 = (4/5 * phys.E_f * fixed.P_e)/(4 * pi^2 * fixed.eta_e * (phys.E_f + phys.E_Li)*sqrt((1 + res.kappa^2) / 2) * fixed.P_w * res.a);

res.B_0 = coil.B_max * (res.R_0 - res.a - res.b) / res.R_0;

xi = coil.B_max^2 / (4 * phys.mu_0 * phys.E_s * coil.e_crit);
sigma = phys.E_s * coil.e_crit;
nu = coil.B_max / (phys.mu_0 * coil.j_max) * (1 - coil.E/phys.E_s);

res.c = (nu + 2 * xi * (res.kappa * res.a + res.b))/(1-xi);

volumeperwatt = sqrt(2 / (1 + res.kappa^2)) * 4/5 * phys.E_f * ((res.kappa * res.a + res.b + res.c)*(res.a + res.b + res.c) - res.kappa * res.a^2)/(2 * fixed.eta_e * (phys.E_f + phys.E_Li) * res.a * fixed.P_w);

% Numerical minimization of the cost to find the cost optimizing value of a, while
% maximizing the usage of magnets

res.aval = double(vpasolve(diff(volumeperwatt,a)==0,[0 5]));


% fill in other parameters
res.cval = double(subs(res.c, a, res.aval));
res.R_0val = double(subs(res.R_0, a, res.aval));
res.B_0val = double(subs(res.B_0, a, res.aval));
end

function [res] = plasmaparameters(res, fixed, phys)

% PLASMAPARAMETERS Calculates p, n, tau_E and I_p of the plasma based on the required
% power output, steady state ignition criterium and ELMy H mode scaling law and the predetermined volume of the plasma

res.V_p =  2 * pi^2 * res.R_0 * res.kappa * res.a^2;

res.p = sqrt(16 * fixed.P_e * (phys.k * res.T)^2 / (fixed.eta_e * (phys.E_f + phys.E_Li) * res.V_p * phys.sigma_v));
res.n = res.p/(2 * phys.k * res.T);

res.tau_e = 3/2 * res.p / (1/4 * res.n^2 * phys.sigma_v * 1/5 * phys.E_f);

P_alpha = 1/4 * res.n^2 * phys.sigma_v * 1/5 * phys.E_f * res.V_p;
res.I_p = ((res.tau_e * (P_alpha * 1e-6)^(0.69)) / (0.145 * phys.H * res.R_0^(1.39) * res.a^(0.58) * res.kappa^(0.78) * (res.n*1e-20)^(0.41) * res.B_0^(0.15) * phys.A^(0.19)))^(1/0.93)*1e6;

end

function [] = plasmalimits(as, res, phys)

% PLASMALIMITS Calculates the normalized ratio of the plasma parameter with prominent
% limits. There is compliance if all limits are < 1
syms a

% External kink safety factor

q_star = 2 * pi * res.a^2 * res.B_0 / (phys.mu_0 * res.R_0 * res.I_p) * ((1 + res.kappa^2) / 2);

q_lim = 2/q_star;

% Troyon beta limit

beta = res.p / (res.B_0^2 / (2 * phys.mu_0));

beta_troyon = 2.8 / 100 * (res.I_p * 1e-6 / (res.a * res.B_0));

troyon_lim = beta/beta_troyon;

% Greenwald density limit

n_greenwald = res.I_p * 1e-6 / (pi * res.a^2) * 1e20;

greenwald_lim = res.n / n_greenwald;

% Required neoclassical bootstrap current fraction (profile assumptions for
% every plasma parameter, see pages 525-526 of FB)

f_b = 4/3 * res.kappa^(1.4) * 2.8/100 * q_star / (res.a / res.R_0)^0.5;

bootstrap_lim = 0.75/f_b;

figure(1)
hold on
plot(as,double(subs(bootstrap_lim,a,as)),'-r','DisplayName','0.75/f_b');
plot(as,subs(greenwald_lim,a,as),'--r','DisplayName','n/n_G');
plot(as,subs(troyon_lim,a,as),'-.r','DisplayName','{\beta}/{\beta_T}');
plot(as,subs(q_lim,a,as),':r','DisplayName','2/q_{*}');  
legend
hold off

end
