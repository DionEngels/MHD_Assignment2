function [c] = CoilForceBalance(a, R_0, phys, coil, res)
%COILFORCEBALANCE Calcualtes the force balance to find c
syms c B_max
F_m = 2 * pi * R_0 * B_max^2 / phys.mu_0 * (2 * a + 2 * res.b + c) / 2;
F_t = coil.sigma_max * 4 * pi * R_0 * c;
c = solve(F_m == F_t, c);
end