% -------------------------------------------------------------------------
% VALIDATION: TP_TEA
%
% Compute equilibrium composition at defined temperature and pressure.
% Reproduce the example case of TEA by Jasmina Blecic.
% URL RESULTS TEA:
% https://github.com/dzesmin/TEA/tree/master/doc/examples/quick_example/results 
%   
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Nov 30 2021
% -------------------------------------------------------------------------
%% INITIALIZE
LS = {'C', 'CH4', 'CO2', 'CO', 'H2', 'H', 'H2O', 'He', 'N2', 'N', 'NH3', 'O'};
self = App(LS);
%% GET INDEXES OF SPECIES
index = find_ind(self.S.LS, LS);
%% TUNNING PARAMETERS
self.TN.tolN = 1e-25;
%% INITIAL CONDITIONS
load Validation_TP_TEA Pressure Temp n_H n_He n_C n_N n_O results_TEA

T = linspace(Temp(1), Temp(end), 500);
p = logspace(-5, 2, 500);

for i=length(T):-1:1
    self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
    self.PD.S_Fuel     = {'H', 'He', 'C', 'N', 'O'};
    self.PD.N_Fuel     = [n_H, n_He, n_C, n_N, n_O];
    %% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'pP', p(i), 'TP', T(i)); 
    %% SOLVE PROBLEM
    self = SolveProblem(self, 'TP');
    %% SAVE RESULTS
    molar_fractions_CT(:, i) = [self.PS.strP{1}.Xi(index(1));...
                        self.PS.strP{1}.Xi(index(2));...
                        self.PS.strP{1}.Xi(index(3));...
                        self.PS.strP{1}.Xi(index(4));...
                        self.PS.strP{1}.Xi(index(5));...
                        self.PS.strP{1}.Xi(index(6));...
                        self.PS.strP{1}.Xi(index(7));...
                        self.PS.strP{1}.Xi(index(8));...
                        self.PS.strP{1}.Xi(index(9));...
                        self.PS.strP{1}.Xi(index(10));...
                        self.PS.strP{1}.Xi(index(11));...
                        self.PS.strP{1}.Xi(index(12))]; 
end
%% GET RESULTS TEA FROM results_TEA
molar_fractions_TEA = [results_TEA.n_C,...
    results_TEA.n_CH4, results_TEA.n_CO2, results_TEA.n_CO, results_TEA.n_H2,...
    results_TEA.n_H, results_TEA.n_H2O, results_TEA.n_He, results_TEA.n_N2,...
    results_TEA.n_N, results_TEA.n_NH3, results_TEA.n_O]';
%% DISPLAY RESULTS - VALIDATION (PLOT)
nfrec = 3;
self.C.mintol_display = self.TN.tolN;
f = figure;
set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
    'outerposition',[0.05 0.05 0.9 0.9]);
axes = gca;
set(axes,'LineWidth',self.Misc.config.linewidth,'FontSize',self.Misc.config.fontsize-2,'BoxStyle','full')
grid(axes, 'off'); box(axes, 'off'); hold(axes, 'on'); axes.Layer = 'Top';
xlabel(axes, 'Temperature, T [K]','FontSize',self.Misc.config.fontsize,'interpreter','latex');
ylabel(axes, 'Molar fraction, $X_i$','FontSize',self.Misc.config.fontsize,'interpreter','latex');
set(axes,'yscale','log')
xlim(axes, [min(Temp), max(Temp)])
ylim(axes, [self.C.mintol_display, 1])
colorbw = brewermap(length(LS), 'Spectral');
for i=1:length(LS)
    plot(T, molar_fractions_CT(i, :), 'LineWidth', self.Misc.config.linewidth, 'color', colorbw(i,:));
end
for i=1:length(LS)
    plot(Temp(1:nfrec:end), molar_fractions_TEA(i, 1:nfrec:end), 'd', 'LineWidth', self.Misc.config.linewidth, 'color', colorbw(i,:));
end
legendname = {'C', 'CH$_4$', 'CO$_2$', 'CO', 'H$_2$', 'H', 'H$_2$O', 'He', 'N$_2$', 'N', 'NH$_3$', 'O',...
              'TEA: C', 'TEA: CH$_4$', 'TEA: CO$_2$', 'TEA: CO', 'TEA: H$_2$', 'TEA: H', 'TEA: H$_2$O',...
              'TEA: He', 'TEA: N$_2$', 'TEA: N', 'TEA: NH$_3$', 'TEA: O'};
legend(legendname, 'FontSize', self.Misc.config.fontsize-6, 'Location', 'northeastoutside', 'interpreter', 'latex');
