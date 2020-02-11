clearvars -except strThProp strR strP
% alpha = [0,5,10,20,50,100];
alpha = [0,5];
gamma = strR{1}.cP/strR{1}.cV;
Gamma = (gamma+1)/(gamma-1);
aux_gamma = 1/Gamma;
v = linspace(-aux_gamma*2,5,300);
aux1 = ones(1,length(v));

for i=length(alpha):-1:1
    p(i,:) = (2*alpha(i)+Gamma-v)./(Gamma.*v-1);
%     mu(i,:) = (p(i,:)-1)./(v-1);
end
%% FIGURES
% Plot configuration
fpath = 'D:\Google Drive\Tesis\Thesis\Figures\';
color = colours;
linewidth = 1.6;
fontsize = 24;
components = ['C','H','O','N'];
f = figure;
set(f,'units','normalized','innerposition',[0 0.1 1 0.9],...
        'outerposition',[0 0.1 1 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
hold on; grid on; box on;
xlim([0,5])
ylim([-aux_gamma*2,5])
xlabel("Relaci\'on vol\'umenes, $v_2/v_1$",'FontSize',fontsize+4,'interpreter','latex');
ylabel("Relaci\'on presiones, $p_2/p_1$",'FontSize',fontsize+4,'interpreter','latex');

plot(v(1,30:end),p(1,30:end),'LineWidth',linewidth,'color',color(2,:))
plot(v(1,30:end),p(2,30:end),'LineWidth',linewidth,'color',color(2,:))
plot(aux1,v,'LineWidth',linewidth-0.4,'color','k')
plot(v,aux1,'LineWidth',linewidth-0.4,'color','k')
plot(aux_gamma.*aux1,v,'--','LineWidth',linewidth,'color','k')
plot(v,-aux_gamma.*aux1,'--','LineWidth',linewidth,'color','k')