% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_REFLECTED
%
% Compute shock polar plots at T1 = 226.65 K and p1 = 0.0117 bar
% (altitude 30 km), a set of 39 species considered, an initial
% shock front Mach numbers = 20, and a deflection angle theta = 35;
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
% Last update Jan 10 2023
% -------------------------------------------------------------------------
    
LS = {'eminus', 'Ar', 'Arplus', 'N', ...
      'Nplus', 'Nminus', 'NO', 'NOplus', 'NO2', ...
      'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
      'N2minus', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
      'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
      'O2minus', 'O3'};

% LS = {'N2', 'O2', 'Ar'};

%% INITIALIZE
self = App(LS);

self.TN.FLAG_FAST = true;
self.TN.FLAG_TCHEM_FROZEN = false;
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 226.65, 'pR', 0.0117);
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar'};
self.PD.N_Oxidizer = [78, 21, 1] ./ 21;
sound_velocity = compute_sound(self.PD.TR.value, self.PD.pR.value, self.PD.S_Oxidizer, self.PD.N_Oxidizer, 'self', self);
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
Mach_number = 20;
self = set_prop(self, 'u1', sound_velocity * Mach_number, 'theta', 35);
%% SOLVE PROBLEM
self = solve_problem(self, 'SHOCK_POLAR_R');
%% DISPLAY RESULTS (PLOTS)
post_results(self);