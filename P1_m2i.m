%% P1 - COMBUSTION M2i: thermochemistry and kinetics
%% Exercise 1: Molar fraction of the reactants 
app.PD.S_Fuel = {'C2H6'};      % Fuel
app.PD.phi.Value = 0.5:0.1:1; % Equivalence ratio
app.PD.CompleteOrIncomplete = 'incomplete'; 
app.PD.TR.Value = 700;         % Temperature of the reactants [K]
app.PD.pR.Value = 1;           % Pressure of the reactants [bar]
p1_ex1_m2i; % Calculations
%% Exercise 2: Molar fractions of the products
displaysweepresults(app.PS.strP,app.PD.phi.Value,app.S.NameSpecies,app.C.mintol_display);
%% Exercise 3: Adiabatic flame temperature
p1_ex3_m2i;
%% Exercise 4: Molar fraction of NO at the output of the chamber
p = 1*1e5;           % pressure [Pa]
X_NO_0 = 0;          % Initial molar fraction of NO
X_O_0  = 5.8115e-06; % Initial molar fraction of O
X_O2_0 = 7.9912e-02; % Initial molar fraction of O2
X_N2_0 = 7.5205e-01; % Initial molar fraction of N2
X_N_0  = 0;          % Initial molar fraction of N
X_OH_0 = 2.2051e-03; % Initial molar fraction of OH
X_H_0  = 0;          % Initial molar fraction of H
% Initial concentration of the chemical species
t0 = 0;              % Initial time [s]
tf = 3e-4;           % Final time   [s]
R = 8.314472;        % Ideal gas constant [J/mol K]
for i=app.C.l_phi:-1:1
    C0 = p/(R*T(i))*[X_NO_0;X_O_0;X_O2_0;X_N2_0;X_N_0;X_OH_0;X_H_0]; 
    [Case{i}.t,Case{i}.C,Case{i}.Xi] = NO_Zeldovich(T(i),p,t0,tf,C0);% [NO,O,O2,N2,N,OH,H]
end
%% FIGURES
% Plot configuration
% fpath = 'D:\Google Drive\Tesis\Slides\AAUgraphics\chemical_equilibrium\';
color = colours;
linewidth = 2.4;
fontsize = 24;
leg = {'NO','O','O2','N2','N','OH','H'};
j=1;

f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
    'outerposition',[0.1 0.1 0.9 0.8])
set(axes,'LineWidth',linewidth-0.2,'FontSize',fontsize+2,'BoxStyle','full')
grid on; box on; hold on
xlim([t(1),t(end)])
ylim([-inf 1])
xlabel('Time, $t^*$','FontSize',fontsize+4,'interpreter','latex');
ylabel("Molar fraction, $X_i$",'FontSize',fontsize+4,'interpreter','latex');

d(1) = plot(t,Xi(:,1),'color',color(1,:),'LineWidth',linewidth-0.4);
d(2) = plot(t,Xi(:,2),'color',color(2,:),'LineWidth',linewidth-0.4);
d(3) = plot(t,Xi(:,3),'color',color(3,:),'LineWidth',linewidth-0.4);
d(4) = plot(t,Xi(:,4),'color',color(4,:),'LineWidth',linewidth-0.4);
d(5) = plot(t,Xi(:,5),'color',color(5,:),'LineWidth',linewidth-0.4);
d(6) = plot(t,Xi(:,6),'color',color(6,:),'LineWidth',linewidth-0.4);
d(7) = plot(t,Xi(:,7),'LineWidth',linewidth-0.4);
legend(leg,'FontSize',fontsize,'Location','northeastoutside','interpreter','latex');
set(gca,'yscale','log')



    