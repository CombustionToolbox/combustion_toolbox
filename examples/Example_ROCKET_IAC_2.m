% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Infinite-Area-Chamber (IAC)
%
% Compute rocket propellant performance considering an Infinite-Area-Chamber
% for lean to rich MHF3-MON25 mixtures at 101.325 bar, a set of 24 species
% considered, a set of equivalence ratios phi contained in (2, 5) [-], and
% an exit area ratio A_exit/A_throat = 3 
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Feb 08 2024
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot Formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 298.15, 'pR', 100 * 1.01325, 'phi', 2);
self.PD.S_Fuel     = {'CH6N2bLb', 'N2H4bLb'};
self.PD.N_Fuel     = convert_weight_percentage_to_moles(self.PD.S_Fuel, [86, 14], self.DB);
self.PD.S_Oxidizer = {'N2O4bLb', 'N2O3'};
self.PD.wt_ratio_oxidizers = [36.67, 63.33];
self.PD.FLAG_IAC = true;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'Aratio', 3);
%% SOLVE PROBLEM
self = solve_problem(self, 'ROCKET');
%% DISPLAY RESULTS (PLOTS)
post_results(self);