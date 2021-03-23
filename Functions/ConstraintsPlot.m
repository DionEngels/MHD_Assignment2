function [] = ConstraintsPlot(limit, res, settings)
%CONSTRAINTSPLOT Plots the regime in which the constraints are valid

if settings.use_rebco
    magnet = 'ReBaCuO';
else
    magnet = 'Nb3Sn';
end

figure(1)
title(strcat("Constraints against ", char(settings.var), " for " ,magnet, " material"))
hold on
plot(settings.values,double(subs(limit.bootstrap,settings.var,settings.values)),'-b','DisplayName','0.75/f_b');
plot(settings.values,subs(limit.greenwald,settings.var,settings.values),'--b','DisplayName','n/n_G');
plot(settings.values,subs(limit.troyon,settings.var,settings.values),'-.b','DisplayName','{\beta}/{\beta_T}');
plot(settings.values,subs(limit.q,settings.var,settings.values),':b','DisplayName','2/q_{*}');  

ylim([0, 3.5])
legend
hold off
ylabel(sprintf("Parameter relative to constraint.\n < 1 is fulfilled"));
switch settings.variable
    case 1
        xlabel(sprintf("%s [m]", char(settings.var)))
    case 2
        xlabel(sprintf("%s [W]", char(settings.var)))
    case 3
        xlabel(sprintf("%s [W]", char(settings.var)))
    case 4
        xlabel(sprintf("%s [-]", char(settings.var)))
end


if settings.plotvolumeperwatt
    figure(2)
    plot(settings.values,double(subs(res.volumeperwatt/(10^-6),settings.var,settings.values)),'-b');
    title(strcat("Volume per Watt for " ,magnet, " material"))
    ylabel(sprintf("Volume per Watt [m^3/MW]"));
    switch settings.variable
        case 1
            xlabel(sprintf("%s [m]", char(settings.var)))
        case 2
            xlabel(sprintf("%s [W]", char(settings.var)))
        case 3
            xlabel(sprintf("%s [W]", char(settings.var)))
        case 4
            xlabel(sprintf("%s [-]", char(settings.var)))
    end
end

end

