function [] = ConstraintsPrint(res, limits, q, f_b)
%CONSTRAINTSPRINT Prints if the result fix the constraints or not
LogicalStr = {'false', 'true'};

constraint_kappa = res.kappa < limits.kappa;
constraint_safety = q > limits.q;
constraint_troyon = res.beta < limits.troyon;
constraint_greenwald = res.n < limits.greenwald;
constraint_bootstrap = f_b > limits.bootstrap;

n_constraints_passed = constraint_kappa + constraint_safety + constraint_troyon + constraint_greenwald + constraint_bootstrap;
fprintf('Kappa: %d < %d. Fulfilled: %s\n', res.kappa, limits.kappa, LogicalStr{constraint_kappa + 1})
fprintf('Safety factor: %d < %d. Fulfilled: %s\n', limits.q, q, LogicalStr{constraint_safety + 1})
fprintf('Troyon limit: %d < %d. Fulfilled: %s\n', res.beta, limits.troyon, LogicalStr{constraint_troyon + 1})
fprintf('Greenwald limit: %d < %d. Fulfilled: %s\n', res.n, limits.greenwald, LogicalStr{constraint_greenwald + 1})
fprintf('Bootstrap current: %d < %d. Fulfilled: %s\n', limits.bootstrap, f_b, LogicalStr{constraint_bootstrap + 1})
fprintf('Result: %u passed, %u failed\n', n_constraints_passed, 5 - n_constraints_passed)
end

