color = colours;
linewidth = 2;
fontsize = 24;

f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.7 0.9],...
    'outerposition',[0.1 0.1 0.7 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
hold on; grid on; grid minor; box on;

xlabel('Equivalence ratio, $\phi$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Temperature, $T [K]$','FontSize',fontsize+10,'interpreter','latex');
plot(app.PD.phi.Value,T,'LineWidth',linewidth,...
        'MarkerFaceColor',color(1,:),'MarkerEdgeColor','black',...
        'MarkerSize',fontsize-12,'color',color(1,:));