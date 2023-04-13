% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_REFLECTED_THERMO
%
% Compute shock polar plots at T1 = 226.65 K and p1 = 0.0117 bar
% (altitude 30 km), a set of 28 species considered, an initial
% shock front Mach number = 20, deflection angle theta = 35, and different
% thermochemical models (chemical equilibrium and frozen chemistry).
%    
% LS== {'eminus', 'Ar', 'Arplus', 'N', ...
%       'Nplus', 'Nminus', 'NO', 'NOplus', 'NO2', ...
%       'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
%       'N2minus', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
%       'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
%       'O2minus', 'O3'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Apr 13 2023
% -------------------------------------------------------------------------
    
LS = {'eminus', 'Ar', 'Arplus', 'N', ...
      'Nplus', 'Nminus', 'NO', 'NOplus', 'NO2', ...
      'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
      'N2minus', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
      'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
      'O2minus', 'O3'};

% Initialize
self = App(LS);
% Miscellaneous
self.TN.FLAG_FAST = true;
self.Misc.FLAG_RESULTS = false;
if FLAG_FROZEN
    self.Misc.config.linestyle = '--';
end
% Thermochemical model
self.PD.FLAG_FROZEN = false;
% Initial conditions
self = set_prop(self, 'TR', 226.65, 'pR', 0.0117);
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar'};
self.PD.N_Oxidizer = [78, 21, 1] ./ 21;
% Additional inputs (depends of the problem selected)
self = set_prop(self, 'M1', 20, 'theta', 35);
% Solve problem
self = solve_problem(self, 'SHOCK_POLAR_R');
% Display results (plots)
post_results(self);