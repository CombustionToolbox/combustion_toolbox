% -------------------------------------------------------------------------
% EXAMPLE: EXOPLANET WASP-43b - METALLICITY 1
%
% Compute equilibrium vertical composition with a metallicity 1 of WASP-43b
%   
% URL RESULTS TEA:
% https://github.com/dzesmin/RRC-BlecicEtal-2015a-ApJS-TEA/tree/master/Fig6/WASP43b-solar
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Oct 12 2022
% -------------------------------------------------------------------------

LS = {'C2H2_acetylene', 'C2H4', 'C', 'CH4', 'CO2', 'CO', 'H2', 'H2O', 'H2S', 'H', 'HCN', 'He', 'HS_M', 'N2', 'N', 'NH3', 'O', 'S'};
%% INITIALIZE
self = App(LS);
%% INITIAL CONDITIONS
metallicity = 1;
Fuel = {'H', 'He', 'C', 'N', 'O', 'S'};
Ni_abundances = abundances2moles(Fuel, 'abundances.txt', metallicity)';
T = linspace(100, 4000, 300);
p = logspace(-5, 2, 300);

self.PD.S_Fuel = Fuel;
self.PD.N_Fuel = Ni_abundances;
self = set_prop(self, 'TR', 300, 'pR', 1);
self = set_prop(self, 'TP', T, 'pP', p);
%% SOLVE PROBLEM
self = solve_problem(self, 'TP');
%% POSTPROCESSING CONFIGURATION
self.Misc.config.label_type = 'long';
%% DISPLAY RESULTS (PLOTS)
plot_molar_fractions(self, self.PS.strP, 'Xi', 'p', 'ydir', 'reverse', 'xscale', 'log');