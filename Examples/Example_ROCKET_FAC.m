% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Infinite-Area-Chamber (IAC)
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (2, 5) [-]
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update April 12 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('HYDROGEN_L');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 90, 'pR', 100 * 1.01325, 'phi', 2);
self.PD.S_Fuel     = {'H2bLb'};
self.PD.S_Oxidizer = {'O2bLb'};
self.PD.FLAG_IAC = false;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'Aratio_c', 2);
%% SOLVE PROBLEM
self = SolveProblem(self, 'ROCKET');
%% DISPLAY RESULTS (PLOTS)
postResults(self);