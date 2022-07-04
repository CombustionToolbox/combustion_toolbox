% -------------------------------------------------------------------------
% EXAMPLE: TP
%
% Compute equilibrium composition of a blast furnace at defined temperature
% 1050 K and pressure 1.01325 bar. A set of 10 gaseous species and 2 
% condensed species are considered
%
% Reference:
% W.D. Madeley & J.M. Toguri (1973) The application of free energy
% minimization techniques to determine equilibrium compositions in systems
% of metallurgical interest, Canadian Metallurgical Quarterly, 12:1, 71-78,
% DOI: 10.1179/cmq.1973.12.1.71
%   
% Species == {'O2', 'N2', 'H2O', 'CH4', 'CO', 'CO2', 'H2', ...
%             'CHO_M', 'CH2O_M', 'OH', 'Febab', 'CaObcrb'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Jul 04 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App({'O2', 'N2', 'H2O', 'CH4', 'CO', 'CO2', 'H2', 'CHO_M', 'CH2O_M', 'OH', 'Febab', 'CaObcrb'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 1050, 'pR', 1 * 1.01325);
self.PD.S_Fuel     = {'O2', 'N2', 'H2O', 'CH4', 'Fe3O4bcrb', 'Febab', 'Cbgrb', 'CaCO3bcrb', 'CaObcrb'};
self.PD.N_Fuel     = [20.46, 187.1, 1.775, 2.2554, 13.1, 3.527, 85.59, 0.1499, 0.6063];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 1050); 
self.TN.tolN = 1e-30;
%% SOLVE PROBLEM
self = SolveProblem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
%% GET RESULTS
moles_CT = moles(self.PS.strP{1});
moles_RAND = [1.793e-19; 1.871e2; 2.377e-1; 1.97e-3; 8.143e1; 6.865; 6.641; 9.46e-7; 1.264e-7; 8.411e-11; 4.283e1; 7.562e-1; 0; 0; 0];
%% DISPLAY TABLE SPECIES AND MOLES
T = table(self.S.LS', moles_CT, moles_RAND, moles_CT - moles_RAND);
T.Properties.VariableNames = {'Species', 'Moles', 'Moles RAND', 'Error'}

