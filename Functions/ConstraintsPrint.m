function [] = ConstraintsPrint(res, limits)
%CONSTRAINTSPRINT Prints if the result fix the constraints or not
LogicalStr = {'false', 'true'};
syms a

lim.q = double(subs(limits.q, a, res.aval));
lim.troyon = double(subs(limits.troyon, a, res.aval));
lim.greenwald = double(subs(limits.greenwald, a, res.aval));
lim.bootstrap = double(subs(limits.bootstrap, a, res.aval));

constraint_safety = lim.q < 1;
constraint_troyon = lim.troyon < 1;
constraint_greenwald = lim.greenwald < 1;
constraint_bootstrap = lim.bootstrap < 1;

n_constraints_passed = constraint_safety + constraint_troyon + constraint_greenwald + constraint_bootstrap;
fprintf('Safety factor: %d. Fulfilled: %s\n', lim.q, LogicalStr{constraint_safety + 1})
fprintf('Troyon limit: %d. Fulfilled: %s\n', lim.troyon,LogicalStr{constraint_troyon + 1})
fprintf('Greenwald limit: %d. Fulfilled: %s\n', lim.greenwald, LogicalStr{constraint_greenwald + 1})
fprintf('Bootstrap current: %d. Fulfilled: %s\n', lim.bootstrap, LogicalStr{constraint_bootstrap + 1})
fprintf('Result: %u passed, %u failed\n', n_constraints_passed, 4 - n_constraints_passed)
end

