function plot_density_YFuel(strR,strR_Fuel)
linewidth = 2;
fontsize = 24;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
        'outerposition',[0.05 0.05 0.9 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
movegui(f,'center')
grid on; box on; hold on; axis tight
xlabel('$Y_{Fuel}\ [-]$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Density, $\rho\ [kg/m^3]$','FontSize',fontsize+10,'interpreter','latex');

Nstruct = length(strR);
for i=Nstruct:-1:1
    YFuel(i) = strR_Fuel.mi/strR{i}.mi;
    rho(i)   = strR{i}.rho;
end
plot(YFuel,rho,'LineWidth',linewidth);