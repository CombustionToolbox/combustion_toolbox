%% VALIDATION: HP SOOT FORMATION 
clc; close all; clear
savefigure = false; % Save figure == true
final = true;
filename = {"JFM_p1","JFM_p10","JFM_p20"};
%% LOAD CASES COMBUSTION-TOOLBOX
TR = 300:50:700;                   % Temperature of reactants
cases = 3;                  % Number of cases to evaluate
NE = length(TR);
colorbw = brewermap(NE,'Spectral');
                  
prop(1) = read_CT('ch4_1_T.mat'); phi_c(1) = 2;
prop(2) = read_CT_strR('ch4_1_T_strR.mat'); phi_c(2) = 2;
% prop(1) = read_CT('h2_1_T.mat'); phi_c(1) = 2;
% prop(2) = read_CT_strR('h2_1_T_strR.mat'); phi_c(2) = 2;
% prop(3) = read_CT('ch4_2_T.mat'); phi_c(1) = 2;
% prop(4) = read_CT_strR('ch4_2_T_strR.mat'); phi_c(2) = 2;

prop(3) = read_CT('ch4_3_T.mat'); phi_c(1) = 2;
prop(4) = read_CT_strR('ch4_3_T_strR.mat'); phi_c(2) = 2;
% prop(3) = read_CT('h2_2_T.mat'); phi_c(1) = 2;
% prop(4) = read_CT_strR('h2_2_T_strR.mat'); phi_c(2) = 2;

prop(5) = read_CT('ch4_4_T.mat'); phi_c(1) = 2;
prop(6) = read_CT_strR('ch4_4_T_strR.mat'); phi_c(2) = 2;
% prop(5) = read_CT('h2_3_T.mat'); phi_c(1) = 2;
% prop(6) = read_CT_strR('h2_3_T_strR.mat'); phi_c(2) = 2;

phi2(1).case = 1;
%% FIGURES
% Plot configuration
% fpath = 'C:\Users\user\Google Drive\Tesis\poster MSC11\Figures\';
color = colours;
linewidth = 2;
fontsize = 25;
language = 2; % 1 == spanish, other number -> english
%% Heat release overdriven
clear q leg
colorbw = brewermap(length(app.PD.phi.Value),'Spectral');

for i=length(app.PD.phi.Value):-1:1
    for j=length(app.PD.overdriven):-1:1
        q(i,j)= app.PS.strP{i}.overdriven{j}.DhT*1e3/app.PS.strP{i}.overdriven{j}.mi-app.PS.strR{i}.overdriven{j}.DhT*1e3/app.PS.strR{i}.overdriven{j}.mi...
            +0.5*(app.PS.strP{i}.overdriven{j}.u^2-app.PS.strR{i}.overdriven{j}.u^2);
    end
end
q= q./q(:,1);
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.8 0.9],...
        'outerposition',[0.05 0.05 0.8 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; hold on; axis tight
ylim([0.4,1])
xlabel('Overdriven $U_o/U_{CJ}$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Heat release, $q^*$','FontSize',fontsize+10,'interpreter','latex');
if final
    set(gca,'color','none')
end
tit ='Overdriven Detonation: CH4 + air (21\%O2, 79\%N2)';
% tit ='Overdriven Detonation: C3H8 + air (21\%O2, 79\%N2)';
% tit ='Overdriven Detonation: C2H5OH + air (21\%O2, 79\%N2)';
title(tit,'Interpreter','latex','FontSize',fontsize+10);


for i=length(app.PD.phi.Value):-1:1
    if app.PD.phi.Value(i)< 1
    dl(i) = plot(app.PD.overdriven,q(i,:),'--','LineWidth',linewidth,'Color',colorbw(i,:));
    else
    dl(i) = plot(app.PD.overdriven,q(i,:),'LineWidth',linewidth,'Color',colorbw(i,:));
    end
    leg(i) = {strcat('$\phi = ',sprintf('%.2f',app.PD.phi.Value(i)),'$')};
%     if savefigure
%         saveas(f,filename{i},'epsc');
%     end
end

legend(dl,leg,'FontSize',fontsize,'Location','northeastoutside','interpreter','latex');
%% Heat release

for i=1:2:cases*2
    if i == 1
        load ZND.mat
    elseif i == 3
        load ZND_10.mat
    elseif i == 5
        load ZND_20.mat
    end
    q(i,:)= prop(i).prop(9,:)./(prop(i).prop(9,:)-1).*prop(i).prop(7,:)*1e5./prop(i).prop(2,:)...
        -prop(i+1).prop(9,:)./(prop(i+1).prop(9,:)-1).*prop(i+1).prop(7,:)*1e5./prop(i+1).prop(2,:)...
        +0.5*(prop(i).prop(8,:).^2-prop(i+1).prop(8,:).^2);
    %     +prop(i).prop(10,:)*1e3-prop(i+1).prop(10,:)*1e3...
    
    q2(i,:)= prop(i).prop(13,:)*1e3-prop(i+1).prop(13,:)*1e3...
        +0.5*(prop(i).prop(8,:).^2-prop(i+1).prop(8,:).^2);
    q_ZND(i,:)= prop_ZND{1,1}.prop(13,:)*1e3-prop_ZND{1,1+1}.prop(13,:)*1e3...
        +0.5*(prop_ZND{1,1}.prop(8,:).^2-prop_ZND{1,1+1}.prop(8,:).^2);
        for k=length(prop_ZND_all{1,1}.struct):-1:1
                q_ZND_all(i).struct{k}.q(:) = prop_ZND_all{1,1}.struct{k}.prop(13,:)*1e3-prop_ZND{1,1+1}.prop(13,k)*1e3*ones(1,length(prop_ZND_all{1,1}.struct{k}.prop(13,:)))...
                    +0.5*(prop_ZND_all{1,1}.struct{k}.prop(8,:).^2-prop_ZND{1,1+1}.prop(8,k).^2.*ones(1,length(prop_ZND_all{1,1}.struct{k}.prop(13,:))));
        end
    for j=length(TR):-1:1
        x_ZND{i,j}.pos = prop_ZND_all{1,1}.struct{j}.prop(14,:);
    end
end
%% MOLE FRACTION
% Plot configuration
nfrec = 3;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.8 0.9],...
        'outerposition',[0.05 0.05 0.8 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; hold on; axis tight
xlabel('Temperature $T\ [K]$','FontSize',fontsize+10,'interpreter','latex');
ylabel('$h_P - h_R\ [kJ/kg]$','FontSize',fontsize+10,'interpreter','latex');
if final
        set(gca,'color','none')
end
tit ='CJ Detonation: CH4 + air (21\%O2, 79\%N2)';
title(tit,'Interpreter','latex','FontSize',fontsize+10);
i=0;
for j=1:2:cases*2 
    i=i+1;
    dl(j) = plot(prop(2).prop(1,:),-(prop(j).prop(10,:).*1e3-prop(j+1).prop(10,:).*1e3)*1e-3,'LineWidth',linewidth,'Color',colorbw(i,:));
end

leg = {'$p_R = 1\ bar$';'$p_R = 10\ bar$';'$p_R = 20\ bar$'};
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');

if savefigure
    saveas(f,filename{j},'epsc');
end
%%
% Plot configuration
nfrec = 3;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.8 0.9],...
        'outerposition',[0.05 0.05 0.8 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; hold on; axis tight
xlabel('Temperature, $T_R\ [K]$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Heat release, $q\ [kJ/kg]$','FontSize',fontsize+10,'interpreter','latex');
if final
        set(gca,'color','none')
end
tit ='CJ Detonation: CH4 + air (21\%O2, 79\%N2)';
title(tit,'Interpreter','latex','FontSize',fontsize+10);
i=0;
for j=1:2:cases*2
    i=i+1;
    dl2(j) = plot(prop(2).prop(1,:),q2(j,:)*1e-3,'LineWidth',linewidth,'Color',colorbw(i,:));
end

leg = {'$p_R = 1\ bar$';'$p_R = 10\ bar$';'$p_R = 20\ bar$'};
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');

if savefigure
    saveas(f,filename{j},'epsc');
end


%%
% Plot configuration
nfrec = 3;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.8 0.9],...
        'outerposition',[0.05 0.05 0.8 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; hold on; axis tight
xlabel('Temperature, $T_R\ [K]$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Heat release, $q\ [kJ/kg]$','FontSize',fontsize+10,'interpreter','latex');
if final
        set(gca,'color','none')
end
tit ='ZND structure: CH4 + air (21\%O2, 79\%N2)';
title(tit,'Interpreter','latex','FontSize',fontsize+10);
i=0;
for j=1:2:cases*2
    i=i+1;
    dl2(j) = plot(prop(2).prop(1,1:end-1),q_ZND(j,1:end-1)*1e-3,'LineWidth',linewidth,'Color',colorbw(i,:));
end

leg = {'$p_R = 1\ bar$';'$p_R = 10\ bar$';'$p_R = 20\ bar$'};
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');

if savefigure
    saveas(f,filename{j},'epsc');
end
%% HEAT REALEASE ZND ALL - SAME PLOT
% mov=VideoWriter('q_ZND_all','MPEG-4');
% set(mov,'FrameRate',1);
% open(mov);
% set(gca,'nextplot','replacechildren');
leg_p{1} = '$p_R = 1\ bar$';
leg_p{3} = '$p_R = 10\ bar$';
leg_p{5} = '$p_R = 20\ bar$';
% t=0;
for j=1:2:cases*2
%     t=t+1;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.85 0.9],...
        'outerposition',[0.05 0.05 0.85 0.9]);
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; hold on; axis tight
xlabel('Distance, $x^*$','FontSize',fontsize+10,'interpreter','latex');
ylabel('Heat release, $q^*$','FontSize',fontsize+10,'interpreter','latex');
xlim([0,1])
ylim([0,1])
if final
        set(gca,'color','none')
end
tit ='ZND structure: CH4 + air (21\%O2, 79\%N2)';
title({tit;leg_p{1,j}},'Interpreter','latex','FontSize',fontsize+10);
for i=1:length(TR)-1
    dl3(i) = plot(x_ZND{j,i}.pos./x_ZND{j,1}.pos(end),q_ZND_all(j).struct{i}.q(:)./q_ZND_all(5).struct{1}.q(end),'LineWidth',linewidth,'Color',colorbw(i,:));
    dl4(i) = plot(x_ZND{j,i}.pos(end)/x_ZND{j,1}.pos(end),q_ZND_all(j).struct{i}.q(end)./q_ZND_all(5).struct{1}.q(end),'d',...
        'LineWidth',linewidth,'MarkerFaceColor',colorbw(i,:),'MarkerEdgeColor','black',...
            'MarkerSize',fontsize-16,'Color',colorbw(i,:));
end

leg = {'$T_R = 300\ K$';'$T_R = 350\ K$';'$T_R = 400\ K$';'$T_R = 450\ K$';...
    '$T_R = 500\ K$';'$T_R = 550\ K$';'$T_R = 600\ K$';'$T_R = 650\ K$';'$T_R = 700\ K$'};
legend(dl3,leg,'FontSize',fontsize,'Location','northwest','interpreter','latex');

if savefigure
    saveas(f,filename{j},'epsc');
end
% if t<4
% MM(t)=getframe(f);
% end
end
% writeVideo(mov,MM);
% close(mov);
% %% HEAT REALEASE ZND ALL - MOVIE
% mov=VideoWriter('q_ZND_all','MPEG-4');
% set(mov,'FrameRate',1);
% open(mov);
% set(gca,'nextplot','replacechildren');
% for i=1:length(TR)
%     fig = figure;
%     set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
%     grid on; box on;hold on; axis tight
%     xlim([0,1e-3])
%     ylim([0,3000])
% %     shading interp; 
%     tit = ['$Temperature = $',num2str(TR(i))];
%     plot(prop_ZND_all{1,1}.struct{i}.prop(14,:),q_ZND_all(1).struct{i}.q(:)*1e-3,'LineWidth',linewidth,'Color',colorbw(1,:))
% %     plot(prop_ZND_all{1,1}.struct{1}.prop(14,:),q_ZND_all(3).struct{k}.q(:)*1e-3,'LineWidth',linewidth,'Color',colorbw(3,:))
% %     plot(prop_ZND_all{1,1}.struct{1}.prop(14,:),q_ZND_all(5).struct{k}.q(:)*1e-3,'LineWidth',linewidth,'Color',colorbw(5,:))
% 
%     title(tit,'interpreter','latex','FontSize',20);
%     set(fig,'Visible','off');
%     MM(i)=getframe(fig);
% %     filename2 = strcat(filename,'-',num2str(i),'.png');
% %     saveas(fig,fullfile(fpath,filename2));
% end
% writeVideo(mov,MM);
% % movie(M);
% close(mov);