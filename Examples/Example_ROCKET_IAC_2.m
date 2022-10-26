% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Infinite-Area-Chamber (IAC)
%
% Compute rocket propellant performance considering an Infinite-Area-Chamber
% for lean to rich MON25-MHF3 mixtures at 101.325 bar, a set of 24 species
% considered, a set of equivalence ratios phi contained in (2, 5) [-], and
% an exit area ratio A_exit/A_throat = 3 
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Jul 30 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot Formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 298.15, 'pR', 100 * 1.01325, 'phi', 2);
self.PD.S_Fuel     = {'N2O4bLb', 'N2O3'};
self.PD.N_Fuel     = convert_weight_percentage_to_moles(self.PD.S_Fuel, [36.67, 63.33], self.DB);
self.PD.S_Oxidizer = {'CH6N2bLb', 'N2H4bLb'};
self.PD.wt_ratio_oxidizers = [86, 14];
self.PD.FLAG_IAC = true;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'Aratio', 3);
%% SOLVE PROBLEM
self = solve_problem(self, 'ROCKET');
%% DISPLAY RESULTS (PLOTS)
post_results(self);