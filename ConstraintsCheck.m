function [] = ConstraintsCheck(res, settings)
%CONSTRAINTSCHECK Function that checks the design to see if they fit the constriants
q = (5 * res.a^2 * res.kappa * res.B_0) / (res.R_0 * res.I_p);
f_b = 1.3 * res.kappa^(1/4) * res.beta * q / res.epsilon^(1/2);

limits.kappa = 2;
limits.q = 2;
limits.troyon = 0.028 * res.I_p / res.a / res.B_0;
limits.greenwald = res.I_p / (pi * res.a^2);
limits.bootstrap = 0.8;

if settings.constraint_plot
    ConstraintsPlot(res, limits, q, f_b, settings);
else
    ConstraintsPrint(res, limits, q, f_b);
end