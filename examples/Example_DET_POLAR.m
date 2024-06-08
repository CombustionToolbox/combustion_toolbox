% -------------------------------------------------------------------------
% EXAMPLE: DET_POLAR
%
% Compute detonation polar curves at standard conditions, a set of 30 species
% considered, and a set of pre-shock Mach numbers M1 = [5, 6, 7, 8, 10], or
% what is the same, drive factors M1/M1_cj = [1.0382, 1.2458, 1.4534,...
% 1.6611, 2.0763]
%    
% Hydrogen == {'H2O','H2','O2','N2','He','Ar','H','HNO',...
%              'HNO3','NH','NH2OH','NO3','N2H2','N2O3','N3','OH',...
%              'HNO2','N','NH3','NO2','N2O','N2H4','N2O5','O','O3',...
%              'HO2','NH2','H2O2','N3H','NH2NO2'}
%
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Oct 07 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Hydrogen');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1);
self.PD.S_Fuel = {'H2'};
self.PD.S_Oxidizer = {'N2', 'O2'};
self.PD.ratio_oxidizers_O2 = [79, 21]/21;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'drive_factor', [1.0382, 1.2458, 1.4534, 1.6611, 2.0763]);
%% TUNING PROPERTIES
self.TN.N_points_polar = 300; % Number of points to compute polar curves
%% SOLVE PROBLEM
self = solve_problem(self, 'DET_POLAR');
%% DISPLAY RESULTS (PLOTS)
post_results(self);