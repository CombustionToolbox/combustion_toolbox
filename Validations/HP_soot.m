%% VALIDATION: HP SOOT FORMATION 
clc; close all; clear
savefigure = false; % Save figure == true
final = true;
filename = {"soot_benzene_T","soot_methane_T","soot_acetylene_T"};
%% LOAD CASES COMBUSTION-TOOLBOX
TR = 300;                   % Temperature of reactants
cases = 3;                  % Number of cases to evaluate
load Data_NameSpecies.mat;  % data from Combustion-Toolbox
% Species to check                        
% species = {'CN','H','Cbgrb','C2H2_acetylene','CO2','CO','C2N2','HNC','HCN','H2','N2','OH','H2O','O2'};
species = {'Cbgrb','CO2','CO','HNC','HCN','H2','OH','H2O','O2','CH4'};
NE = length(species);
colorbw = brewermap(NE,'Spectral');
for t=numel(species(1,:)):-1:1
    for k = 1:length(NameSpecies)
        if strcmp(NameSpecies{k},species{t})
            Nt(t) = k;
        end
    end
end
                   
% prop(1) = read_CT('Data_soot_benzene_correct_11_06_2019_strP.mat'); phi_c(1) = 2.5;
prop(1) = read_CT('Data_benzene_msc11.mat'); phi_c(1) = 2.5;
% prop(2) = read_CT('Data_soot_methane_correct_2_strP.mat'); phi_c(2) = 4;
% prop(3) = read_CT('Data_soot_acetylene_correct_30_05_2019_strP.mat'); phi_c(3) = 2.5;
% prop(3) = read_CT('Data_soot_acetylene_correct_11_06_2019_strP.mat'); phi_c(3) = 2.5;
prop(3) = read_CT('Data_acetylene_msc11.mat'); phi_c(3) = 2.5;
% prop(3) = read_CT('Data_soot_acetylene_correct_2_strP');   phi_c(3) = 2.5;
% prop(3) = read_CT('Data_soot_acetylene_zoom_strP'); phi_c(3) = 2.5;
% prop(3) = read_CT('Data_soot_acetylene_zoom_NEW_31_05_2019_strP.mat'); phi_c(3) = 2.5;
% prop(3) = read_CT('Data_soot_acetylene_zoom_NEW_6_06_2019_strP.mat'); phi_c(3) = 2.5;
% prop(3) = read_CT('Data_soot_acetylene_zoom_31_05_2019_strP.mat'); phi_c(3) = 2.5;
%% PROPERTIES NASA
% [Prop(1),s(1),phi_nasa(1).case] = data_CEA('benzene_soot.out',species);
[Prop(1),s(1),phi_nasa(1).case] = data_CEA('benzene_msc11_all.out',species);
% [Prop(2),s(2),phi_nasa(2).case] = data_CEA('methane_soot.out',species);
% [Prop(3),s(3),phi_nasa(3).case] = data_CEA('test_soot_acetylene_zoom.txt',species);
% for i=length(prop(3).X(1,:)):-1:1
%     idx_Xi(:,i) = prop(3).X(:,i)>1e-15;
% end
% idx_Xi = prop(3).X>1e-17;
% for i=length(prop(3).X(:,1)):-1:1
%     idx(i,1) = any(idx_Xi(i,:));
% end
% species=NameSpecies(idx);
% [Prop(3),s(3),phi_nasa(3).case] = data_CEA('test_soot_acetylene.txt',species);
[Prop(3),s(3),phi_nasa(3).case] = data_CEA('acetylene_msc11_all.out',species);
%% PROPERTIES CANTERA + SDT TOOLBOX
mech = 'soot.cti';
% mech = 'gri30_highT.cti';
gas = Solution(mech);
% phi2(1).case = [2:0.01:4];
phi2(1).case = [0.5:0.01:4.5];
for j=length(phi2(1).case):-1:1
    pR = 1*1e5;
    reaction = strcat('c6h6:1,O2:',num2str(7.5/phi2(1).case(j)),',N2:',num2str(7.5*(79/21)/phi2(1).case(j)));
    set(gas,'T',TR,'P',pR,'X',reaction);
    equilibrate(gas,'HP');
    RHO(j) = density(gas);
%     H(j) = enthalpy_mass(gas)*1e-3;
%     E(j) = intEnergy_mass(gas)*1e-3;
    S(j) = entropy_mass(gas)*1e-3;
    T(j) = temperature(gas);
    P(j) = pressure(gas)*1e-5;
    G(j) = gibbs_mass(gas)*1e-3;
end
Prop_cc(1).case = [T;RHO;S;G];
%%
% clear T RHO S G
gas = GRI30;
phi2(2).case = 0.5:0.01:4.5;
% for j=length(phi2(2).case):-1:1
%     pR = 1*1e5;
%     reaction = strcat('CH4:1,O2:',num2str(2/phi2(2).case(j)),',N2:',num2str(2*(79/21)/phi2(2).case(j)));
%     set(gas,'T',TR,'P',pR,'X',reaction);
%     equilibrate(gas,'HP');
%     RHO(j) = density(gas);
% %     H(j) = enthalpy_mass(gas)*1e-3;
% %     E(j) = intEnergy_mass(gas)*1e-3;
%     S(j) = entropy_mass(gas)*1e-3;
%     T(j) = temperature(gas);
% %     P(j) = pressure(gas)*1e-5;
%     G(j) = gibbs_mass(gas)*1e-3;
% end
% Prop_cc(2).case = [T;RHO;S;G];
%%
clear T RHO S G
gas = GRI30;
% phi2(3).case = 2:0.01:3.5;
phi2(3).case = [0.5:0.01:4.5];
% phi2(3).case = [2.5:0.001:2.64]; % zoom
for j=length(phi2(3).case):-1:1
    pR = 1*1e5;
    reaction = strcat('C2H2:1,O2:',num2str(2.5/phi2(3).case(j)),',N2:',num2str(2.5*(79/21)/phi2(3).case(j)));
    set(gas,'T',TR,'P',pR,'X',reaction);
    equilibrate(gas,'HP');
    RHO(j) = density(gas);
%     H(j) = enthalpy_mass(gas)*1e-3;
%     E(j) = intEnergy_mass(gas)*1e-3;
    S(j) = entropy_mass(gas)*1e-3;
    T(j) = temperature(gas);
    P(j) = pressure(gas)*1e-5;
    G(j) = gibbs_mass(gas)*1e-3;
end
Prop_cc(3).case = [T;RHO;S;G;P];
%% FIGURES
% Plot configuration
fpath = 'C:\Users\user\Google Drive\Tesis\poster MSC11\Figures\';
color = colours;
linewidth = 2;
fontsize = 25;
language = 2; % 1 == spanish, other number -> english
%% Plot Properties
% clear leg
% name_save = ["T","rho","s"]; 
% if language == 1
%     Name_prop = ["Densidad, $\rho$","Entalp\'ia, $h$","Energ\'ia interna, $e$","Entrop\'ia, $s$",...
%     "Temperatura, $T$","Presi\'on, $p$","Velocidad de detonaci\'on, $u_1$","Velocidad incidente, $u_2$"];
% else
%     Name_prop = ["Temperature, $T$ ","Density, $\rho$","Entropy, $s$"];
% end
% unit_prop = ["$K$","$kg/m^3$","$kJ/kg-K$"];
% 
% for j=3:3
% for i=3:-1:1
%     f = figure;
%     set(f,'units','normalized','innerposition',[0.1 0.1 0.7 0.9],...
%         'outerposition',[0.1 0.1 0.7 0.9]);
%     set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
%     hold on; grid on; grid minor; box on;
%     %     xlim([0.5 1.9])
%     if language == 1
%         xlabel('Temperatura, $T_R [K]$','FontSize',fontsize+10,'interpreter','latex');
%     else
%         xlabel('Equivalence ratio, $\phi$','FontSize',fontsize+10,'interpreter','latex');
%     end
%     ylabel(strcat(Name_prop(i),' [',unit_prop(i),']'),'FontSize',fontsize+10,'interpreter','latex');
%     
%     plot(phi2(j).case,prop(j).prop(i,:),'LineWidth',linewidth,...
%         'MarkerFaceColor',color(1,:),'MarkerEdgeColor','black',...
%         'MarkerSize',fontsize-12,'color',color(1,:))
%     plot(phi2(j).case,Prop_cc(j).case(i,:),'LineWidth',linewidth,...
%         'MarkerFaceColor',color(4,:),'MarkerEdgeColor','black',...
%         'MarkerSize',fontsize-12,'color',color(4,:))
%     plot(phi_nasa(j).case,Prop(j).case(i,:),'d','LineWidth',linewidth,...
%         'MarkerFaceColor',color(2,:),'MarkerEdgeColor','black',...
%         'MarkerSize',fontsize-12,'color',color(2,:))
% 
%     %%%% Error
% %     for j=1:1
% %     if j==1  % Last case didnt converge for phi=0.5
% %         error2 = norm(prop(j).prop(i,:)-Prop(j).case(i,1:end-1))/norm(prop(j).prop(i,:))*100;
% %         error3 = norm(prop(j).prop(i,:)-Prop_g(j).case(i,1:end-1))/norm(prop(j).prop(i,:))*100;
% %         error4 = norm(prop(j).prop(i,:)-Prop_cc(j).case(i,1:end-1))/norm(prop(j).prop(i,:))*100;
% %     else
% %         error2 = norm(prop(1).prop(i,:)-Prop(1).case(i,:))/norm(prop(1).prop(i,:))*100;
% %         error3 = norm(prop(1).prop(i,:)-Prop_g(1).case(i,:))/norm(prop(1).prop(i,:))*100;
% %         error4 = norm(prop(1).prop(i,:)-Prop_cc(1).case(i,:))/norm(prop(1).prop(i,:))*100;
% %     end
%     
% %     end
% %     maxerror2 = max(error2);
% %     maxerror3 = max(error3);
% %     maxerror4 = max(error4);
%     %%%%
% %     tit1 = strcat('$\epsilon_{CEA} = ',num2str(maxerror2,'%.2f'),' \%$');
% %     title({tit1},'Interpreter','latex','FontSize',fontsize+4);
%     
% %     tit1 = strcat('$\epsilon_{CEA} = ',num2str(maxerror2,'%.2f'),' \%$',...
% %         '; $\ \epsilon_{CANTERA} = ',num2str(maxerror4,'%.2f'),' \%$;');
% %     tit2 = strcat('$\ \epsilon_{GASEQ} = ',num2str(maxerror3,'%.2f'),'$');
% 
%     tit1(1).tit ='Case HP: Benzene $(C_6H_6)$';
%     tit1(2).tit ='Case HP: Methane $(CH_4)$';
%     tit1(3).tit ='Case HP: Acetylene $(C_2H_2)$';
%     title({tit1(j).tit},'Interpreter','latex','FontSize',fontsize+10);
%     
% %     leg = {'Combustion-Toolbox','CANTERA + SDT','CEA'};
%     leg = {'Combustion-Toolbox','CEA'};
%     legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');
%     %         legend(leg,'FontSize',fontsize,'Location','northoutside','interpreter','latex','NumColumns',4);
%     
%     if savefigure
% %         filename2 = strcat(fpath,strcat('HP_STAMS_',name_save(i)));
% %         filename2 = strcat('legend_slides_',num2str(j),'_',name_save(i));
%         saveas(f,filename{j},'epsc');
%     end
% end
% end

%% MOLE FRACTION
% Plot configuration
% color = colours;
nfrec = 3;
for j=3:cases
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.8 0.9],...
        'outerposition',[0.05 0.05 0.8 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; hold on; axis tight
% xlim([phi_nasa(3).case(1),phi_nasa(3).case(end)])
% ylim([mintol,1])
xlabel('Equivalence Ratio $\phi$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Molar fraction $X_i$','FontSize',fontsize+10,'interpreter','latex');
leg =NameSpecies(Nt);
% cmap = linspecer(numel(species));
cmap = colormap(hsv(numel(species)));
aux_y = logspace(-16,1);
% for j=1:cases
if final
        set(gca,'color','none')
end
l_phic = plot(1*ones(1,length(aux_y)),aux_y,'LineWidth',linewidth,'Color','k');
l_phic = plot(phi_c(j)*ones(1,length(aux_y)),aux_y,'--','LineWidth',linewidth,'Color','k');
for i=1:numel(Nt)
    dl(i) = plot(phi2(j).case,prop(j).X(Nt(i),:),'LineWidth',linewidth,'Color',colorbw(i,:));
    plot(phi_nasa(j).case,s(j).nasa(:,i),'d','LineWidth',linewidth,...
            'MarkerFaceColor',dl(i).Color,'MarkerEdgeColor','black',...
            'MarkerSize',fontsize-16,'color',dl(i).Color)
    if mod(i,nfrec)==0
        loc_label = 'right';
    else
        loc_label = 'left';
    end
    if strcmpi(NameSpecies{Nt(i)},'Cbgrb')
        label(dl(i),'Cgr','FontSize',fontsize,'location',loc_label,'Interpreter','latex');
    else
        label(dl(i),NameSpecies{Nt(i)},'FontSize',fontsize,'location',loc_label,'Interpreter','latex');
    end
end
annotation(f,'textbox',...
    [0.35594652406417 0.35680751173709 0.0431176470588236 0.0481220657276995],...
    'String','CH4',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FitBoxToText','off','color',colorbw(10,:));

annotation(f,'textbox',...
    [0.565171122994652 0.636150234741784 0.0431176470588236 0.0481220657276995],...
    'String',{'Cgr'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FitBoxToText','off','color',colorbw(1,:));
text(2.54,1940,'$\phi_{gr}$','Interpreter','latex','fontsize',fontsize+3) 


set(gca,'yscale','log')
ylim([1e-16 1])
legend(dl,leg,'FontSize',16,'Location','northoutside','interpreter','latex','NumColumns',NE);
if savefigure
    %         filename2 = strcat(fpath,strcat('HP_STAMS_',name_save(i)));
    %         filename2 = strcat('legend_slides_',num2str(j),'_',name_save(i));
    saveas(f,filename{j},'epsc');
end
end
% end


%% ALL T IN ONE PLOT
clear leg
name_save = ["T","rho","s"]; 
if language == 1
    Name_prop = ["Densidad, $\rho$","Entalp\'ia, $h$","Energ\'ia interna, $e$","Entrop\'ia, $s$",...
    "Temperatura, $T$","Presi\'on, $p$","Velocidad de detonaci\'on, $u_1$","Velocidad incidente, $u_2$"];
else
    Name_prop = ["Temperature, $T$ ","Density, $\rho$","Entropy, $s$"];
end
unit_prop = ["$K$","$kg/m^3$","$kJ/kg-K$"];

for j=1:3:3
for i=1:-1:1
    f = figure;
    set(f,'units','normalized','innerposition',[0.1 0.1 0.8 0.7],...
        'outerposition',[0.1 0.1 0.8 0.7]);
    set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
    if final
        set(gca,'color','none')
    end
    hold on; grid on;
    %     xlim([0.5 1.9])
    ylim([1300 2800])
    yticks([1400:200:2800])
    if language == 1
        xlabel('Temperatura, $T_R [K]$','FontSize',fontsize+10,'interpreter','latex');
    else
        xlabel('Equivalence ratio, $\phi$','FontSize',fontsize+10,'interpreter','latex');
    end
    ylabel(strcat(Name_prop(i),' [',unit_prop(i),']'),'FontSize',fontsize+10,'interpreter','latex');
    
    if j==1
        aux_y = linspace(1300,2800);
        l_phic = plot(1*ones(1,length(aux_y)),aux_y,'LineWidth',linewidth,'Color','k');
    end
    if j==1
%         aux_y = linspace(1300,2600);
        aux_y = linspace(1300,2800);
%         aux_y2 = linspace(2740,2800);
        l_phic = plot(phi_c(j)*ones(1,length(aux_y)),aux_y,'--','LineWidth',linewidth,'Color','k');
%         l_phic = plot(phi_c(j)*ones(1,length(aux_y)),aux_y2,'--','LineWidth',linewidth,'Color','k');
        text(2.54,1940,'$\phi_{gr}$','Interpreter','latex','fontsize',fontsize+3) 
        
%         aux_y = linspace(1300,2600);
%         aux_y = linspace(1300,2800);
%         aux_y2 = linspace(2740,2800);
%         l_phic = plot(phi_c(3)*ones(1,length(aux_y)),aux_y,'--','LineWidth',linewidth,'Color','k');
%         l_phic = plot(phi_c(3)*ones(1,length(aux_y)),aux_y2,'--','LineWidth',linewidth,'Color','k');
%         text(2.54,1940,'$\phi_{gr,C2H4}$','Interpreter','latex','fontsize',fontsize+3) 
    end
    plot(phi2(j).case,prop(j).prop(i,:),'LineWidth',linewidth,...
        'MarkerFaceColor',color(1,:),'MarkerEdgeColor','black',...
        'MarkerSize',fontsize-12,'color',color(1,:))
%     plot(phi2(j).case,Prop_cc(j).case(i,:),'LineWidth',linewidth,...
%         'MarkerFaceColor',color(4,:),'MarkerEdgeColor','black',...
%         'MarkerSize',fontsize-12,'color',color(4,:))
    plot(phi_nasa(j).case,Prop(j).case(i,:),'d','LineWidth',linewidth,...
        'MarkerFaceColor',color(2,:),'MarkerEdgeColor','black',...
        'MarkerSize',fontsize-12,'color',color(2,:))
    
    plot(phi2(j+2).case,prop(j+2).prop(i,:),'LineWidth',linewidth,...
        'MarkerFaceColor',color(1,:),'MarkerEdgeColor','black',...
        'MarkerSize',fontsize-12,'color',color(1,:))
%     plot(phi2(j).case,Prop_cc(j).case(i,:),'LineWidth',linewidth,...
%         'MarkerFaceColor',color(4,:),'MarkerEdgeColor','black',...
%         'MarkerSize',fontsize-12,'color',color(4,:))
    plot(phi_nasa(j+2).case,Prop(j+2).case(i,:),'d','LineWidth',linewidth,...
        'MarkerFaceColor',color(2,:),'MarkerEdgeColor','black',...
        'MarkerSize',fontsize-12,'color',color(2,:))
    
    
    %%%% Error
%     for j=1:1
%     if j==1  % Last case didnt converge for phi=0.5
%         error2 = norm(prop(j).prop(i,:)-Prop(j).case(i,1:end-1))/norm(prop(j).prop(i,:))*100;
%         error3 = norm(prop(j).prop(i,:)-Prop_g(j).case(i,1:end-1))/norm(prop(j).prop(i,:))*100;
%         error4 = norm(prop(j).prop(i,:)-Prop_cc(j).case(i,1:end-1))/norm(prop(j).prop(i,:))*100;
%     else
%         error2 = norm(prop(1).prop(i,:)-Prop(1).case(i,:))/norm(prop(1).prop(i,:))*100;
%         error3 = norm(prop(1).prop(i,:)-Prop_g(1).case(i,:))/norm(prop(1).prop(i,:))*100;
%         error4 = norm(prop(1).prop(i,:)-Prop_cc(1).case(i,:))/norm(prop(1).prop(i,:))*100;
%     end
    
%     end
%     maxerror2 = max(error2);
%     maxerror3 = max(error3);
%     maxerror4 = max(error4);
    %%%%
%     tit1 = strcat('$\epsilon_{CEA} = ',num2str(maxerror2,'%.2f'),' \%$');
%     title({tit1},'Interpreter','latex','FontSize',fontsize+4);
    
%     tit1 = strcat('$\epsilon_{CEA} = ',num2str(maxerror2,'%.2f'),' \%$',...
%         '; $\ \epsilon_{CANTERA} = ',num2str(maxerror4,'%.2f'),' \%$;');
%     tit2 = strcat('$\ \epsilon_{GASEQ} = ',num2str(maxerror3,'%.2f'),'$');

%     tit1(1).tit ='Case HP: Benzene $(C_6H_6)$';
%     tit1(2).tit ='Case HP: Methane $(CH_4)$';
%     tit1(3).tit ='Case HP: Acetylene $(C_2H_2)$';
%     title({tit1(j).tit},'Interpreter','latex','FontSize',fontsize+10);
    
%     leg = {'Combustion-Toolbox','CANTERA + SDT','CEA'};
%     leg = {'Combustion-Toolbox','CEA'};
%     legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');
    %         legend(leg,'FontSize',fontsize,'Location','northoutside','interpreter','latex','NumColumns',4);
    
    if savefigure
%         filename2 = strcat(fpath,strcat('HP_STAMS_',name_save(i)));
%         filename2 = strcat('legend_slides_',num2str(j),'_',name_save(i));
        saveas(f,filename{j},'epsc');
    end
end
end