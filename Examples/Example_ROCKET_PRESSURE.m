% -------------------------------------------------------------------------
% EXAMPLE: ROCKET PRESSURE
%
% Compute adiabatic temperature and equilibrium composition for a range of
% pressures (4.83,  34.47, 137.9 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (0.5, 2) [-]
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
% Last update March 14 2022
% -------------------------------------------------------------------------

%% INITIALIZE
% profile on
pressure = [4.83, 34.47, 137.9]; % [bar]
self = App('HYDROGEN_L');
for i = length(pressure):-1:1
    self = App('fast', self.DB_master, self.DB, self.S.LS);
    %% INITIAL CONDITIONS
    self = set_prop(self, 'TR', 90, 'pR', pressure(i), 'phi', 0.5:0.1:2);
    self.PD.S_Fuel     = {'H2bLb'};
    self.PD.S_Oxidizer = {'O2bLb'};
    %% SOLVE PROBLEM
    self = SolveProblem(self, 'ROCKET');
    self_all{i} = self;
end
%% DISPLAY RESULTS (PLOTS)
postResults(self_all);
% profile viewer