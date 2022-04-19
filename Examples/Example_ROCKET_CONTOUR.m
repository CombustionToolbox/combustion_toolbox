% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering a Finite-Area-Chamber (FAC)
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (1, 5) [-]
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


%% SOLVE PROBLEM
phi = 1:0.1:3;
Aratio = 1.1:0.1:2.62;
for i = length(phi):-1:1
    %% INITIALIZE
    self = App('SOOT FORMATION');
    %% INITIAL CONDITIONS
    self = set_prop(self, 'TR', 298.15, 'pR', 22, 'phi', 1.44922);
    self.PD.S_Fuel     = {'RP_1'};
    self.PD.S_Oxidizer = {'O2bLb'};
    self.PD.FLAG_IAC = false;
    %% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'Aratio_c', 2, 'Aratio', Aratio);
    self = set_prop(self, 'TR', 298.15, 'pR', 22, 'phi', phi(i));
    self = SolveProblem(self, 'ROCKET');
    Z(:, i) = cell2vector(self.PS.strP, 'I_sp');
end
%% DISPLAY RESULTS (PLOTS)
% postResults(self);

[X, Y] = meshgrid(phi, Aratio);
contourf(X, Y, Z);