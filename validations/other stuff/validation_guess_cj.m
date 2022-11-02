N = length(app.PS.strP);
for j=N:-1:1
    T_guess(j) = app.PS.strP{1,j}.T_guess;
    p_guess(j) = app.PS.strP{1,j}.p_guess * 1e-5;
    T(j) = app.PS.strP{1,j}.T;
    p(j) = app.PS.strP{1,j}.p;
end

f1 = figure;
set(f1,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
        'outerposition',[0.1 0.1 0.9 0.8])
ax = axes(f1);
hold(ax, 'on');
set(ax,'LineWidth', app.Misc.config.linewidth,'FontSize', app.Misc.config.fontsize-2,'BoxStyle','full')
xlabel(ax, 'Equivalence ratio', 'FontSize', app.Misc.config.fontsize,'interpreter','latex');
ylabel(ax, 'Temperature [K]', 'FontSize', app.Misc.config.fontsize,'interpreter','latex');
plot(app.PD.phi.value, T, 'linewidth', app.Misc.config.linewidth)
plot(app.PD.phi.value, T_guess, 'linewidth', app.Misc.config.linewidth)
plot(app.PD.phi.value, 6*app.PS.strR{1,1}.T*ones(N,1), 'linewidth', app.Misc.config.linewidth)
legend(ax, {'numerical', 'guess', 'guess-old'}, 'FontSize', app.Misc.config.fontsize-2,'Location','northeast','interpreter','latex');

f2 = figure;
set(f2,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
        'outerposition',[0.1 0.1 0.9 0.8])
ax = axes(f2);
hold(ax, 'on');
set(ax,'LineWidth', app.Misc.config.linewidth,'FontSize', app.Misc.config.fontsize-2,'BoxStyle','full')
xlabel(ax, 'Equivalence ratio','FontSize', app.Misc.config.fontsize,'interpreter','latex');
ylabel(ax, 'Pressure [bar]','FontSize', app.Misc.config.fontsize,'interpreter','latex');
plot(app.PD.phi.value, p, 'linewidth', app.Misc.config.linewidth)
plot(app.PD.phi.value, p_guess, 'linewidth', app.Misc.config.linewidth)
plot(app.PD.phi.value, 15*app.PS.strR{1,1}.p*ones(N,1), 'linewidth', app.Misc.config.linewidth)
legend(ax, {'numerical', 'guess', 'guess-old'}, 'FontSize', app.Misc.config.fontsize-2,'Location','northeast','interpreter','latex');