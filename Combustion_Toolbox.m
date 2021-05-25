%% COMBUSTION TOOLBOX
% TYPE OF PROBLEMS:
%
% * TP ------> Equilibrium composition at defined T and p
% * HP ------> Adiabatic T and composition at constant p
% * SP ------> Isentropic compression/expansion to a specified p
% * TV ------> Equilibrium composition at defined T and constanxt v
% * EV ------> Adiabatic T and composition at constant v
% * SV ------> Isentropic compression/expansion to a specified v
% * SHOCK_I -> Planar incident shock wave
% * SHOCK_R -> Planar reflectet shock wave
% * DET -----> Chapman-Jouget Detonation (CJ upper state)
% * DET_OVERDRIVEN -----> Overdriven Detonation
%
% Authors:
%
% * Alberto Cuadra Lara, Universidad Carlos III de Madrid (UC3M)
% * Marcos Vera Coello,  Universidad Carlos III de Madrid (UC3M)
%
% Last update: 25-May-2021 10:27
%% LOAD DATABASES, GLOBAL PARAMETERS AND MINORS PRODUCTS
addpath(genpath(pwd));
% app = App();
app = App('Soot formation');
% app = App('HC/02/N2 EXTENDED');
% app = App('Hydrogen');
% app = App('NASA ALL');
% app = App('Cbgrb'); 
app.Misc.save_Excel = false;
%% REACTION: COMPLETE OR INCOMPLETE
% app.PD.CompleteOrIncomplete = 'complete';
app.PD.CompleteOrIncomplete = 'incomplete'; 
app.TN.factor_c = 1; % factor_c = 1 (default).
%% CHECK SPECIES MINOR PRODUCTS --> SELECTED DATABASE (strThProp)
[app.strThProp,app.E,app.S,app.M,app.C] = check_database(app.strMaster,app.strThProp,app.E,app.S,app.M,app.C);
%% PROBLEM CONDITIONS
app.PD.TR.value = 300;
% app.PD.TR.vector.value = 300:50:700;
app.PD.pR.value = 1;
app.PD.phi.value = 0.1:0.01:2;
%% INITIALIZATION
app = Initialize(app);
%% PROBLEM TYPE
switch app.PD.ProblemType
    case 'TP' % * TP: Equilibrium composition at defined T and p
        app.PD.ProblemType = 'TP';
%         app.PD.TP.value = [300:10:2000];
        app.PD.TP.value = 2000;
    case 'HP' % * HP: Adiabatic T and composition at constant p
        app.PD.ProblemType = 'HP';
        % app.PD.pR.value = logspace(0,2,20); app.PD.phi.value = 1*ones(1,length(app.PD.pR.value));
    case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
        app.PD.ProblemType = 'SP';
        app.PD.pP.value = 10:10:50; app.PD.phi.value = 1*ones(1,length(app.PD.pP.value));
%         app.PD.pP.value = 10*ones(1,length(app.PD.phi.value));
    case 'TV' % * TV: Equilibrium composition at defined T and constant v
        app.PD.ProblemType = 'TV';
        app.PD.TP.value = 4000;
    case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
        app.PD.ProblemType = 'EV';
        app.PD.pR.value = app.PD.pR.value;
        % app.PD.pR.value = logspace(0,2,20); app.PD.phi.value = 1*ones(1,length(app.PD.pR.value));
    case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
        app.PD.ProblemType = 'SV';
        % REMARK! vP_vR > 1 --> expansion, vP_vR < 1 --> compression
        app.PD.vP_vR.value = 0.5:0.1:2; app.PD.phi.value = 1*ones(1,length(app.PD.vP_vR.value));
    case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
        app.PD.ProblemType = 'SHOCK_I';
        app.PD.u1_vector.value = 400:200:2000; app.PD.phi.value = ones(1,length(app.PD.u1_vector.value));
    case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
        app.PD.ProblemType = 'SHOCK_R';
        app.PD.u1_vector.value = 400:50:3000; app.PD.phi.value = ones(1,length(app.PD.u1_vector.value));
    case 'DET' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
        app.PD.ProblemType = 'DET';
%         app.PD.TR_vector.value = app.PD.TR.value;
    case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION
        app.PD.ProblemType = 'DET_OVERDRIVEN';
        app.PD.overdriven.value  = 1:0.1:1.5;
end
%% CONSTANT
app.C.l_phi = length(app.PD.phi.value);
tic
for i=app.C.l_phi:-1:1 % Evading preallocate struct
%% DEFINE FUEL
% app.PD.R_Fuel = 0; app.PD.phi_t = 1; app.PD.Fuel.x = 0; app.PD.Fuel.eps = 1e-1; app.C.FLAG_Fuel = 0;
% app.PD.S_Fuel = {'CH4','C2H6','C3H8'}; app.PD.N_Fuel = [0.85;0.1;0.05]; 
% app.PD.S_Fuel = {'CH4'}; app.PD.N_Fuel = 1;
app.PD.S_Fuel = {'C2H2_acetylene'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C3H8'}; app.PD.N_Fuel = 1;    
% app.PD.S_Fuel = {'C6H6'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C2H4'}; app.PD.N_Fuel = 1;

app = Define_F(app);
%% DEFINE OXIDIZER
app.PD.S_Oxidizer = {'O2'}; app.PD.N_Oxidizer = app.PD.phi_t/app.PD.phi.value(i);
app = Define_O(app);
%% DEFINE DILUENTS/INERTS
app.PD.S_Inert = {'N2'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.value(i)*79/21;
% app.PD.S_Inert = {'N2','Ar'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.value(i).*[78.09/20.95;0.93/20.95];
app = Define_I(app);
%% COMPUTE PROPERTIES
app = Define_FOI(app, i);
%% PROBLEM TYPE
app = SolveProblem(app, i);
%% DISPLAY RESULTS
if ~strcmp(app.PD.ProblemType,'SHOCK_R') && ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
    displayresults(app.PS.strR{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.namespecies);
elseif ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
    displayresults(app.PS.strR{i},app.PS.str2{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.namespecies); % Display all results SHOCK_R
end
end
disp('TIME:')
toc
%% DISPLAY MOLAR FRACTION VS EQUIVALENCE RATIO
app.M.display_species = {};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','N2','O2'};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb'};
closing(app,app.PS.strP,app.PD.phi.value,app.M.display_species,app.Misc.timer_0,app.S.namespecies,app.C.mintol_display,app.PD.ProblemType);
%% EXCEL I/O
ExportExcel(app);