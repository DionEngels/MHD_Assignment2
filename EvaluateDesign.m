function [res] = EvaluateDesign(res, coil, c, R_0, settings)
%EVALUATEDESIGN Evaluates the design for certain coils and results for a
syms a B_max
res.B_max = coil.B_max;
res.c = double(subs(subs(c, a, res.a), B_max, res.B_max));

res.B_0 = ((res.a + res.b + res.c/2) / 2) * res.B_max / res.R_0;
if ~settings.constraint_plot
    res.R_0 = double(subs(R_0, a, res.a));
end
end

