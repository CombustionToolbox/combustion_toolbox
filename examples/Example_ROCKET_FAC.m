% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Finite-Area-Chamber (FAC)
%
% Compute rocket propellant performance considering an Finite-Area-Chamber
% with an area ratio of the combustion chamber A_chamber/A_throat = 2
% for lean to rich LH2-LOX mixtures at 101.325 bar, a set of 24 species
% considered, a set of equivalence ratios phi contained in (2, 5) [-], and
% an exit area ratio A_exit/A_throat = 3 
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Jul 30 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('HYDROGEN_L');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 90, 'pR', 100 * 1.01325, 'phi', 1:0.05:5);
self.PD.S_Fuel     = {'H2bLb'};
self.PD.S_Oxidizer = {'O2bLb'};
self.PD.FLAG_IAC = false;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'Aratio_c', 2, 'Aratio', 3);
%% SOLVE PROBLEM
self = solve_problem(self, 'ROCKET');
%% DISPLAY RESULTS (PLOTS)
post_results(self);