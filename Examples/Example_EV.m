%{
  EXAMPLE: EV

  Compute equilibrium composition at adiabatic temperature (e.g., 3000 K)
  and constant volume for lean to rich CH4-air mixtures at standard
  conditions, a set of 24 species considered and a set of equivalence
  ratios (phi) contained in (0.5, 5) [-]
    
  Soot formation == {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
                     'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                     'NO','HCO','NH2','NH','N','CH','Cbgrb'}
    
  See wiki or ListSpecies() for more predefined sets of species

  @author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Universidad Carlos III de Madrid
                  
  Last update Oct 22 2021
%}
%% INITIALIZE
self = App('Soot formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 0.5:0.01:5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
% No additional data required. The internal energy and volume are constants.
%% SOLVE PROBLEM
self = SolveProblem(self, 'EV');
%% DISPLAY RESULTS (PLOTS)
closing(self);