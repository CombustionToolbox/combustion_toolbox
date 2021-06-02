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
%% INITIALIZATION
addpath(genpath(pwd));
% app = App();
% app = App('Soot formation');
% app = App('HC/02/N2 extended');
% app = App('HC/02/N2 rich');
app = App('Hydrogen');
% app = App('Nasa all');
% app = App('Cbgrb'); 
% app = App('air'); 
% app = App({'Mgbcrb','MgObcrb','NO','O'}); 
%% REACTION: COMPLETE OR INCOMPLETE
% app.PD.CompleteOrIncomplete = 'complete';
app.PD.CompleteOrIncomplete = 'incomplete';
%% TUNNING AND MISCELANEOUS PARAMETERS
app.Misc.save_Excel = false;
app.TN.factor_c = 1; % factor_c = 1 (default).
%% SOLVER: SEGREGATED, GIBBS OR GIBBS_SOOT
app.PD.solver = 'SEGREGATED';
% app.PD.solver = 'GIBBS';
% app.PD.solver = 'GIBBS REDUCED';
% app.PD.solver = 'GIBBS SOOT';
%% CHECK SPECIES MINOR PRODUCTS --> SELECTED DATABASE (strThProp)
[app.strThProp,app.E,app.S,app.M,app.C] = check_database(app.strMaster,app.strThProp,app.E,app.S,app.M,app.C);
%% PROBLEM CONDITIONS
app.PD.TR.value = 300;
% app.PD.TR.vector.value = 300:50:700;
app.PD.pR.value = 1.01325;
app.PD.phi.value = 0.1:0.01:3;
% app.PD.phi.value = 3;
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
        app.PD.pP.value = 10:1:50; app.PD.phi.value = 1*ones(1,length(app.PD.pP.value));
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
        app.PD.vP_vR.value = 0.5:0.01:2; app.PD.phi.value = 1*ones(1,length(app.PD.vP_vR.value));
    case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
        app.PD.ProblemType = 'SHOCK_I';
        app.PD.u1.value = 400:200:2000; app.PD.phi.value = ones(1,length(app.PD.u1.value));
    case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
        app.PD.ProblemType = 'SHOCK_R';
        app.PD.u1.value = 400:50:3000; app.PD.phi.value = ones(1,length(app.PD.u1.value));
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
% app.PD.S_Fuel = {'CH4','C2H6','C3H8'}; app.PD.N_Fuel = [0.85;0.1;0.05]; 
% app.PD.S_Fuel = {'CH4'}; app.PD.N_Fuel = 1;
% app.PD.S_Fuel = {'C2H2_acetylene'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C3H8'}; app.PD.N_Fuel = 1;    
% app.PD.S_Fuel = {'C6H6'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C2H4'}; app.PD.N_Fuel = 1;
app.PD.S_Fuel = {'H2'}; app.PD.N_Fuel = 1;
% app.PD.S_Fuel = {'Mgbcrb'}; app.PD.N_Fuel = 1;
app = Define_F(app);
%% DEFINE OXIDIZER
app.PD.S_Oxidizer = {'O2'}; app.PD.N_Oxidizer = app.PD.phi_t/app.PD.phi.value(i);
% app.PD.S_Oxidizer = {'O2'}; app.PD.N_Oxidizer = 0.5;
app = Define_O(app);
%% DEFINE DILUENTS/INERTS
app.PD.proportion_N2_O2 = 79/21;
app.PD.S_Inert = {'N2'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.value(i) * app.PD.proportion_N2_O2;
% app.PD.S_Inert = {'N2'}; app.PD.N_Inert = app.PD.N_Oxidizer * app.PD.proportion_N2_O2;
app = Define_I(app);
%% COMPUTE PROPERTIES
app = Define_FOI(app, i);
%% PROBLEM TYPE
app = SolveProblem(app, i);
%% DISPLAY RESULTS COMMAND WINDOW
results(app, i);
end
%% DISPLAY RESULTS (PLOTS)
% app.M.display_species = {};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','N2','O2'};
app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb'};
% app.M.display_species = {'Cbgrb','CO2','CO','HNC','HCN','H2','OH','H2O','O2','CH4'};
closing(app);
%% EXCEL I/O
ExportExcel(app);