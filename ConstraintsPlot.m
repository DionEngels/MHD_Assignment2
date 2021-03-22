function [] = ConstraintsPlot(res, limits, q, f_b, settings)
%CONSTRAINTSPLOT Plots the regime in which the constraints are valid
figure;
hold on
fplot(f_b / limits.bootstrap, settings.a_bounds);
fplot(res.n / limits.greenwald, settings.a_bounds);
fplot(res.beta / limits.troyon, settings.a_bounds);
fplot(q / limits.q, settings.a_bounds);
ylim([0, 2])
end

