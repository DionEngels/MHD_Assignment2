close all

%%%% IMPORTANT: convert every parameter and formula to SI units before
%%%% using it in this program %%%%%%

%% Settings

as = linspace(1,1.7,100); % Range of a values to plot and the resolution
settings_constraint_plot = 1;   % 1 for show plots, 0 for calculate check constraints

%% Physical constants

global mu_0 e k

mu_0 = 4*pi*10^-7;     % H/m
e = 1.60217649*10^-19; % C
k = 1.3806504*10^-23;  % J*C/K

%% Fixed parameters

global P_e P_w sigma_v eta_e A H E_f E_Li E_s 

P_e = 1000*10^6;        % W power provided to the grid
P_w = 4*10^6;           % W/m2 max neutron wall loading
sigma_v = 3*10^(-22);   % m3/s (at 15 keV)
eta_e = 0.4;            % - thermal conversion efficiency
A = 2.5;                % - average atomic mass for D-T 
H = 1;                  % - H mode enhancement factor (tau_E is multiplied with this factor)
E_f = 17.6*10^6*e;      % J energy per fusion reaction
E_Li = 4.8*10^6*e;      % J energy per neutron produced via Li-6 breeding in the blanket
E_s = 210*10^9;         % Pa Youngs modulus of the structural material (steel) 

%% Materials data from the superconducting magnets

B_max_Nb = 13;           % T
B_max_Re = 19;           % T
j_max_Nb = 1000e6;       % A/m^2
j_max_Re = 2500e6;       % A/m^2
e_crit_Nb =  0.2/100;    % - (not % !)
e_crit_Re =  0.62/100;   % - (not % !)
E_Nb =  150*10^9;        % Pa
E_Re =  172*10^9;        % Pa


%% Predetermined parameters that aren't optimized

global b kappa T

kappa = 1.7;            % - (Stabilization limit of the n=0 vertical displacement mode)
b = 1.2;                % m (blanket thickness based on nuclear engineering constraints)
T = 15*10^3*e/k;        % K (Temperature at which the fusion power is optimized for a fixed pressure)

%% Parameter to vary

syms a;                 % m

%% Volume and cost of the device from the wall loading and power output contraints (engineering)

[volperwatt_Nb,c_Nb,sigma_struct_Nb,R_0_Nb] = dimensions(a, B_max_Nb, e_crit_Nb, j_max_Nb, E_Nb);
[volperwatt_Re,c_Re,sigma_struct_Re,R_0_Re] = dimensions(a, B_max_Re, e_crit_Re, j_max_Re, E_Re);

% Numerical minimization of the cost to find the cost optimizing value of a, while
% maximizing the usage of magnets

a_opt_Nb = double(vpasolve(diff(volperwatt_Nb,a)==0,[0 5]));
a_opt_Re = double(vpasolve(diff(volperwatt_Re,a)==0,[0 5]));
c_opt_Nb = double(subs(c_Nb,a_opt_Nb));
c_opt_Re = double(subs(c_Re,a_opt_Re));

%% Resulting plasma parameters based on scaling laws, steady state ignition and the required power output

[p_Nb,n_Nb,tau_E_Nb,I_p_Nb] = plasmaparameters(a,R_0_Nb,B_max_Nb);
[p_Re,n_Re,tau_E_Re,I_p_Re] = plasmaparameters(a,R_0_Re,B_max_Re);

%% Checking if the required plasma parameters are feasible through the normalized limits (plasma physics)

[q_lim_Nb,troyon_lim_Nb,greenwald_lim_Nb,bootstrap_lim_Nb] = plasmalimits(a,R_0_Nb,B_max_Nb,I_p_Nb,p_Nb,n_Nb);
[q_lim_Re,troyon_lim_Re,greenwald_lim_Re,bootstrap_lim_Re] = plasmalimits(a,R_0_Re,B_max_Re,I_p_Re,p_Re,n_Re);

%% Plot constraints or calculate constraints

if settings_constraint_plot
    
    figure(1)
    hold on
    plot(as,double(subs(bootstrap_lim_Nb,a,as)),'-r','DisplayName','0.75/f_b Nb3Sn');
    plot(as,subs(greenwald_lim_Nb,a,as),'--r','DisplayName','n/n_G Nb3Sn');
    plot(as,subs(troyon_lim_Nb,a,as),'-.r','DisplayName','{\beta}/{\beta_T} Nb3Sn');
    plot(as,subs(q_lim_Nb,a,as),':r','DisplayName','2/q_{*} Nb3Sn');  
    plot(as,subs(bootstrap_lim_Re,a,as),'-b','DisplayName','0.75/f_b ReBCuO');
    plot(as,subs(greenwald_lim_Re,a,as),'--b','DisplayName','n/n_G ReBCuO');
    plot(as,subs(troyon_lim_Re,a,as),'-.b','DisplayName','{\beta}/{\beta_T} ReBCuO');
    plot(as,subs(q_lim_Re,a,as),':b','DisplayName','2/q_{*} ReBCuO');
    legend
    hold off
    
    figure(2)
    as = linspace(1,3,100); % Range of a values to plot and the resolution
    hold on
    plot(as,subs(volperwatt_Nb/(10^(-6)),a,as),'-r','DisplayName','Volume per Watt Nb3Sn [m^3/(MW)]');
    plot(as,subs(volperwatt_Re/(10^(-6)),a,as),'-b','DisplayName','Volume per Watt ReBCuO [m^3/(MW)]');
    legend
    hold off
    
else
    % To be reexamined !
    
    LogicalStr = {'false', 'true'};
   
    constraint_safety_Nb = q_lim_Nb < 1;
    constraint_safety_Re = q_lim_Re < 1;
    constraint_troyon_Nb = troyon_lim_Nb < 1;
    constraint_troyon_Re = troyon_lim_Re < 1;
    constraint_greenwald_Nb = greenwald_lim_Nb < 1;
    constraint_greenwald_Nb = greenwald_lim_Re < 1;
    constraint_bootstrap = f_b > limit_bootstrap;
    
    n_constraints_passed = constraint_kappa + constraint_safety + constraint_troyon + constraint_greenwald + constraint_bootstrap;
    fprintf('Kappa constraint is fulfilled: %s\n', LogicalStr{constraint_kappa + 1})
    fprintf('Safety factor constraint is fulfilled: %s\n', LogicalStr{constraint_safety + 1})
    fprintf('Troyon limit constraint is fulfilled: %s\n', LogicalStr{constraint_troyon + 1})
    fprintf('Greenwald limit constraint is fulfilled: %s\n', LogicalStr{constraint_greenwald + 1})
    fprintf('Bootstrap current constraint is fulfilled: %s\n', LogicalStr{constraint_bootstrap + 1})
    fprintf('Result: %u passed, %u failed\n', n_constraints_passed, 5 - n_constraints_passed)
end

function [volumeperwatt,c,sigma,R_0] = dimensions(a, B_max, e_crit, j_max, E_sc)
% COSTANDVOLUME Calculates the mayor radius, plasma volume, total cost, the total magnet thickness and the normal
% stress to which the structural material in the TF field coils is subjected
syms cost c R_0 V_p

global mu_0 b kappa P_e P_w eta_e E_f E_Li E_s 

R_0 = (4/5*E_f*P_e)/(4*pi^2*eta_e*(E_f+E_Li)*sqrt((1+kappa^2)/2)*P_w*a);

B_0 = B_max*(R_0-a-b)/R_0;

xi = B_max^2/(4*mu_0*E_s*e_crit);
sigma = E_s*e_crit;
nu = B_max/(mu_0*j_max)*(1-E_sc/E_s);

c = (nu+2*xi*(kappa*a+b))/(1-xi);

volumeperwatt = sqrt(2/(1+kappa^2))*4/5*E_f*((kappa*a+b+c)*(a+b+c)-kappa*a^2)/(2*eta_e*(E_f+E_Li)*a*P_w);

end

function [p,n,tau_E,I_p] = plasmaparameters(a,R_0,B_max)

% PLASMAPARAMETERS Calculates p, n, tau_E and I_p of the plasma based on the required
% power output, steady state ignition criterium and ELMy H mode scaling law and the predetermined volume of the plasma

global k P_e sigma_v eta_e A H E_f E_Li kappa b T 

V_p =  2*pi^2*R_0*kappa*a^2;

p = sqrt(16*P_e*(k*T)^2/(eta_e*(E_f+E_Li)*V_p*sigma_v));
n = p/(2*k*T);

B_0 = B_max*(R_0-a-b)/R_0;

tau_E = 3/2*p/(1/4*n^2*sigma_v*1/5*E_f);

P_alpha = 1/4*n^2*sigma_v*1/5*E_f*V_p;
I_p = ((tau_E*(P_alpha*10^(-6))^(0.69))/...
    (0.145*H*R_0^(1.39)*a^(0.58)*kappa^(0.78)*(n*10^(-20))^(0.41)*B_0^(0.15)*A^(0.19)))^(1/0.93)*10^6;

end

function [q_lim,troyon_lim,greenwald_lim,bootstrap_lim] = plasmalimits(a,R_0,B_max,I_p,p,n)

% PLASMALIMITS Calculates the normalized ratio of the plasma parameter with prominent
% limits. There is compliance if all limits are < 1

global mu_0 b kappa

% External kink safety factor

B_0 = B_max*(R_0-a-b)/R_0;

q_star = 2*pi*a^2*B_0/(mu_0*R_0*I_p)*((1+kappa^2)/2);

q_lim = 2/q_star;

% Troyon beta limit

beta = p/(B_0^2/(2*mu_0));

beta_troyon = 2.8/100*(I_p*10^(-6)/(a*B_0));

troyon_lim = beta/beta_troyon;

% Greenwald density limit

n_greenwald = I_p*(10^-6)/(pi*a^2)*(10^20);

greenwald_lim = n/n_greenwald;

% Required neoclassical bootstrap current fraction (profile assumptions for
% every plasma parameter, see pages 525-526 of FB)

f_b = 4/3*kappa^(1.4)*2.8/100*q_star/(a/R_0)^0.5;

bootstrap_lim = 0.75/f_b;

end
