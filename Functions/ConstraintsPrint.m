function [] = ConstraintsPrint(limits)
%CONSTRAINTSPRINT Prints if the result fix the constraints or not
LogicalStr = {'false', 'true'};

constraint_safety = limits.q < 1;
constraint_troyon = limits.troyon < 1;
constraint_greenwald = limits.greenwald < 1;
constraint_bootstrap = limits.bootstrap < 1;

n_constraints_passed = constraint_safety + constraint_troyon + constraint_greenwald + constraint_bootstrap;
fprintf('Safety factor: %d. Fulfilled: %s\n', limits.q, LogicalStr{constraint_safety + 1})
fprintf('Troyon limit: %d. Fulfilled: %s\n', limits.troyon,LogicalStr{constraint_troyon + 1})
fprintf('Greenwald limit: %d. Fulfilled: %s\n', limits.greenwald, LogicalStr{constraint_greenwald + 1})
fprintf('Bootstrap current: %d. Fulfilled: %s\n', limits.bootstrap, LogicalStr{constraint_bootstrap + 1})
fprintf('Result: %u passed, %u failed\n', n_constraints_passed, 4 - n_constraints_passed)
end

