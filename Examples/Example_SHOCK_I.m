% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_I
%
% Compute pre-shock and post-shock state for a planar incident shock wave
% at standard conditions, a set of 20 species considered and a set of
% initial shock front velocities (u1) contained in (360, 20000) [m/s]
%    
% Air == {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4',...
%         'N3','C','CO','CO2','Ar','H2O','H2','H','He'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Oct 22 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Air');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 0.5:0.01:5);
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
u1 = logspace(2, 5, 500); u1 = u1(u1<20000); u1 = u1(u1>=360);
self = set_prop(self, 'u1', u1, 'phi', self.PD.phi.value(1) * ones(1, length(u1)));
%% SOLVE PROBLEM
self = SolveProblem(self, 'SHOCK_I');
%% DISPLAY RESULTS (PLOTS)
closing(self);