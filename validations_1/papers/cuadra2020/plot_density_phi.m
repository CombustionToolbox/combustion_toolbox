function [rho, YFuel, W, YFuel_st, phi] = plot_density_phi(strR,strR_Fuel)
linewidth = 2;
fontsize = 24;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
        'outerposition',[0.05 0.05 0.9 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
movegui(f,'center')
grid on; box on; hold on; axis tight
xlabel('$\phi\ [-]$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Density, $\rho\ [kg/m^3]$','FontSize',fontsize+10,'interpreter','latex');

Nstruct = length(strR);
for i=Nstruct:-1:1
    phi(i) = strR{i}.phi;
    YFuel(i) = strR_Fuel.mi/strR{i}.mi;
    rho(i) = strR{i}.rho;
    if strR{i}.phi == 1
        YFuel_st = YFuel(i);
    end
end
if ~exist('YFuel_st')
   YFuel_st = 0; 
end
plot(phi,rho,'LineWidth',linewidth);
W = plot_WorH_YFuel(strR,strR_Fuel,YFuel,rho,true,phi);