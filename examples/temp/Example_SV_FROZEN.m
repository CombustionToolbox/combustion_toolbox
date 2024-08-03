% -------------------------------------------------------------------------
% EXAMPLE: SV FROZEN
% Compute Isentropic compression/expansion and equilibrium composition at 
% a defined set of volume ratios (0.5, 2) for a lean CH4-air mixture at
% 700 K and 10 bar, frozen chemistry, and a equivalence ratio phi 0.5 [-]
%   
% LS == {'CH4', 'O2', 'N2', 'Ar', 'CO2'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update July 22 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App({'CH4', 'O2', 'N2', 'Ar', 'CO2'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 700, 'pR', 10, 'phi', 0.5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.ratio_oxidizers_O2 = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'vP_vR', 0.5:0.01:2); 
%% SOLVE PROBLEM
self = solve_problem(self, 'SV');
%% DISPLAY RESULTS (PLOTS)
post_results(self);