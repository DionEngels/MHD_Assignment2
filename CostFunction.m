function [a] = CostFunction(a, b, c, R_0, fixed, coil, settings)
%COSTFUNCTION Determines the optimum value of a using the cost function
% Also plots
syms B_max

cost_function = 2 * pi^2 * R_0 * ((a + b + c)^2 - a^2) / (fixed.P_e/1e6);
cost_function = subs(cost_function, B_max, coil.B_max);

fplot(cost_function, settings.a_bounds);
ylim([0, 4]);

a = fminbnd(matlabFunction(cost_function), settings.a_bounds(1), settings.a_bounds(2));

end

