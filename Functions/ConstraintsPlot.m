function [] = ConstraintsPlot(limit, settings)
%CONSTRAINTSPLOT Plots the regime in which the constraints are valid

figure;
title(sprintf("Constraints against %s", char(settings.var)))
hold on
plot(settings.values,double(subs(limit.bootstrap,settings.var,settings.values)),'-r','DisplayName','0.75/f_b');
plot(settings.values,subs(limit.greenwald,settings.var,settings.values),'--r','DisplayName','n/n_G');
plot(settings.values,subs(limit.troyon,settings.var,settings.values),'-.r','DisplayName','{\beta}/{\beta_T}');
plot(settings.values,subs(limit.q,settings.var,settings.values),':r','DisplayName','2/q_{*}');  

ylim([0, 4])
legend
hold off
box on
ylabel(sprintf("Parameter relative to constraint.\n >1 is good"));
switch settings.variable
    case 1
        xlabel(sprintf("%s [m]", char(settings.var)))
    case 2
        xlabel(sprintf("%s [W]", char(settings.var)))
end

end

