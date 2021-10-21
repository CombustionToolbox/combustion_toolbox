%{ 
COMBUSTION TOOLBOX @v0.3.7

Type of problems:
    * TP -----------------> Equilibrium composition at defined T and p
    * HP -----------------> Adiabatic T and composition at constant p
    * SP -----------------> Isentropic compression/expansion to a specified p
    * TV -----------------> Equilibrium composition at defined T and constant v
    * EV -----------------> Adiabatic T and composition at constant v
    * SV -----------------> Isentropic compression/expansion to a specified v
    * SHOCK_I ------------> Planar incident shock wave
    * SHOCK_R ------------> Planar reflected shock wave
    * DET ----------------> Chapman-Jouget Detonation (CJ upper state)
    * DET_OVERDRIVEN -----> Overdriven Detonation    
    
    
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Universidad Carlos III de Madrid
                  
Last update Oct 11 2021
---------------------------------------------------------------------- 
%}
addpath(genpath(pwd));

%% INITIALIZE
self = App('Soot formation');
% self = App('HC/02/N2');
% self = App('HC/02/N2 extended');
% self = App('HC/02/N2 rich');
% self = App('HC/02/N2 propellants');
% self = App('Ideal_air');
% self = App('Hydrogen_l');
% self = App({'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4','N3', ...
%     'eminus', 'Nminus', 'Nplus', 'NOplus', 'NO2minus', 'NO3minus', 'N2plus', 'N2minus', 'N2Oplus', ...
%      'Oplus', 'Ominus', 'O2plus', 'O2minus'});
%% PROBLEM CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 0.5:0.01:5);
% * DEFINE FUEL | default: self.PD.N_Fuel = 1
self.PD.S_Fuel = {'CH4'};
%  self.PD.N_Fuel = 1;
% * DEFINE OXIDIZER | default: self.PD.N_Oxidizer = fun(phi == equivalence ratio)
self.PD.S_Oxidizer = {'O2'};
% self.PD.N_Oxidizer = 2;
% * DEFINE DILUENTS/INERTS | default: self.PD.N_Inert = fun(phi == equivalence ratio)
self.PD.S_Inert = {'N2'};
% self.PD.N_Inert = 7.52;
%% PROBLEM TYPE
switch self.PD.ProblemType
    case 'TP' % * TP: Equilibrium composition at defined T and p
        self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 3000);
    case 'HP' % * HP: Adiabatic T and composition at constant p
        self = set_prop(self, 'pP', self.PD.pR.value);
    case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
        pP = 10:1:50;
        self = set_prop(self, 'pP', pP, 'phi', self.PD.phi.value(1) * ones(1, length(pP)));
%         self = set_prop(app, 'pP', 10*ones(1, length(self.PD.phi.value)));
    case 'TV' % * TV: Equilibrium composition at defined T and constant v
        self = set_prop(self, 'TP', 2000, 'pP', self.PD.pR.value); % Pressure value is a guess!
    case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
        self = set_prop(self, 'pP', self.PD.pR.value);
    case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
        % REMARK! vP_vR > 1 --> expansion, vP_vR < 1 --> compression
        vP_vR = 0.5:0.01:2;
        self = set_prop(self, 'vP_vR', vP_vR, 'phi', self.PD.phi.value(1) * ones(1, length(vP_vR)));
    case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
        u1 = logspace(2, 5, 500);
        u1 = u1(u1<20000); u1 = u1(u1>=360);
%         u1 = [356,433,534,658,811,1000,1233,1520,1874,2310,2848,3511,4329,5337,6579,8111,9500,12328,15999,18421,21210,24421,28118,32375,37276,42919,49417,56899,65513];
%         u1 = linspace(360, 9000, 1000);
%         u1 = 20000;
        self = set_prop(self, 'u1', u1, 'phi', self.PD.phi.value(1) * ones(1, length(u1)));
    case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
        u1 = linspace(400, 6000, 1000);
        self = set_prop(self, 'u1', u1, 'phi', self.PD.phi.value(1) * ones(1, length(u1)));
    case 'DET' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
%         self.PD.TR_vector.value = self.PD.TR.value;
    case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION
        overdriven = 1:0.1:10;
        self = set_prop(self, 'overdriven', overdriven, 'phi', self.PD.phi.value(1) * ones(1, length(overdriven)));
end
%% SOLVE SELECTED PROBLEM
tic
self = SolveProblem(self);
toc
%% DISPLAY RESULTS (PLOTS)
% Uncomment to specify a custom set of species to be displayed in the molar fraction figure
% self.Misc.display_species = {'CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb'};
closing(self);