function [] = ConstraintsPlot(res, limits, q, f_b, settings)
%CONSTRAINTSPLOT Plots the regime in which the constraints are valid
figure;
title("Constraints against minor radius")
hold on
fplot(f_b / limits.bootstrap, settings.a_bounds, 'DisplayName', 'f_{b} / f_{limit}');
fplot(res.n / limits.greenwald, settings.a_bounds, 'DisplayName', 'n/n_G');
fplot(res.beta / limits.troyon, settings.a_bounds, 'DisplayName', '{\beta}/{\beta_T}');
fplot(q / limits.q, settings.a_bounds, 'DisplayName', 'q_{*} / q_{limit}');
ylim([0, 2])
legend
end

