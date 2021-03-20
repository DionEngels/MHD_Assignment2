function [c] = CoilForceBalance(R_0, B_max, mu_0, a, b, sigma_max)
%COILFORCEBALANCE Calcualtes the force balance to find c
syms c
F_m = 2 * pi * R_0 * B_max^2 / mu_0 * (2 * a + 2 * b + c) / 2;
F_t = sigma_max * 4 * pi * R_0 * c;
c = solve(F_m == F_t, c);
end

