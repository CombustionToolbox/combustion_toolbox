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
% Last update: 05-Feb-2020 18:07
%% LOAD DATABASES AND GLOBAL PARAMETERS
addpath(genpath(pwd));
app = Initialize(); 
app.Misc.save_Excel = false;
%% REACTION: COMPLETE OR INCOMPLETE
% Specify the type of the reaction: complete or incomplete (dissociation)

% app.PD.CompleteOrIncomplete = 'complete';
app.PD.CompleteOrIncomplete = 'incomplete'; 
app.TN.factor_c = 1; % factor_c = 1 (default).
%% MINORS PRODUCTS
% Specify the minority products to be considered in the product mixture (P) in
% addition to the major species (CO2, CO, H2O, H2, O2, N2, C(gr))
% Also, He and Ar are always included in the database.

app = MinorProducts(app, 'Soot formation');
% app = MinorProducts(app, 'HC/02/N2 EXTENDED');
% app = MinorProducts(app, 'Hydrogen');
% app = MinorProducts(app, 'NASA ALL');
% app = MinorProducts(app, 'Cbgrb');
%% CHECK SPECIES MINOR PRODUCTS --> SELECTED DATABASE (strThProp)
% Checks if any specie of the app.M.minor_products considered are not included
% in the selected database strThProp. In case some is missing compute it.
% IF THE MINOR-SPECIE IS NOT DECLARE AS REACTANT THE THEORY WILL FAIL!!
[app.strThProp,app.E,app.S,app.M,app.C] = check_database(app.strMaster,app.strThProp,app.E,app.S,app.M,app.C);
%% PROBLEM CONDITIONS
% Specify the pressure (app.PD.pR.Value [bar]) and the equivalence ratio (app.PD.phi.Value) of the
% reactant mixture (R). The latter may not be required for non-reacting
% systems, e.g., air dissociation.
% Also, depending the transformation considered it is necessary to specify
% additional parameters.
%
% EXAMPLE:
% * app.PD.TR.Value = 300;         % temperature of the reactive species [K]
% * app.PD.pR.Value = 1;           % pressure of the reactive species [bar]
% * app.PD.phi.Value = 0.4:0.05:2; % equivalence ratio [-]
app.PD.TR.Value = 300;
% app.PD.TR.vector.Value = 300:50:700;
% app.PD.TR.Value = 300;
app.PD.pR.Value = 100;
% app.PD.phi.Value = 2*ones(1,171);
% app.PD.phi.Value = 3:0.1:4;
app.PD.phi.Value = 0.1:0.01:2.4;
% app.PD.phi.Value = 1;
%% INITIALIZATION
[app.E,app.S,app.M,app.C,app.Misc,Problem_selected] = Initialize_2(app.E,app.S,app.M,app.C,app.Misc);
%% PROBLEM TYPE
switch Problem_selected
    case 'TP' % * TP: Equilibrium composition at defined T and p
        app.PD.ProblemType = 'TP';
%         app.PD.TP_vector.Value = [300:10:2000];
        app.PD.TP_vector.Value = 2000;
    case 'HP' % * HP: Adiabatic T and composition at constant p
        app.PD.ProblemType = 'HP';
        app.PD.pR_vector.Value = app.PD.pR.Value;
        % app.PD.pR_vector.Value = logspace(0,2,20); app.PD.phi.Value = 1*ones(1,length(app.PD.pR_vector.Value));
    case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
        app.PD.ProblemType = 'SP';
        app.PD.pP_vector.Value = 10:10:50; app.PD.phi.Value = 1*ones(1,length(app.PD.pP_vector.Value));
%         app.PD.pP_vector.Value = 10*ones(1,length(app.PD.phi.Value));
    case 'TV' % * TV: Equilibrium composition at defined T and constant v
        app.PD.ProblemType = 'TV';
        app.PD.TP_vector.Value = 4000;
    case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
        app.PD.ProblemType = 'EV';
        app.PD.pR_vector.Value = app.PD.pR.Value;
        % app.PD.pR_vector.Value = logspace(0,2,20); app.PD.phi.Value = 1*ones(1,length(app.PD.pR_vector.Value));
    case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
        app.PD.ProblemType = 'SV';
        % REMARK! vP_vR > 1 --> expansion, vP_vR < 1 --> compression
        app.PD.vP_vR_vector.Value = 0.5:0.1:2; app.PD.phi.Value = 1*ones(1,length(app.PD.vP_vR_vector.Value));
    case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
        app.PD.ProblemType = 'SHOCK_I';
        app.PD.u1_vector.Value = 400:200:2000; app.PD.phi.Value = ones(1,length(app.PD.u1_vector.Value));
    case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
        app.PD.ProblemType = 'SHOCK_R';
        app.PD.u1_vector.Value = 400:50:3000; app.PD.phi.Value = ones(1,length(app.PD.u1_vector.Value));
    case 'DET' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
        app.PD.ProblemType = 'DET';
%         app.PD.TR_vector.Value = app.PD.TR.Value;
    case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION
        app.PD.ProblemType = 'DET_OVERDRIVEN';
        app.PD.overdriven  = 1:0.1:1.5;
end
%% CONSTANT
app.C.l_phi = length(app.PD.phi.Value);
tic
for i=app.C.l_phi:-1:1 % Evading preallocate struct
% waitbar(1-i/app.C.l_phi,f,strcat('Case ',sprintf(' %d',app.C.l_phi-i+1),' -',sprintf(' %d',app.C.l_phi)));
% if getappdata(f,'canceling'), break, end
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
app.PD.S_Oxidizer = {'O2'}; app.PD.N_Oxidizer = app.PD.phi_t/app.PD.phi.Value(i);
app = Define_O(app);
%% DEFINE DILUENTS/INERTS
app.PD.S_Inert = {'N2'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.Value(i)*79/21;
% app.PD.S_Inert = {'N2','Ar'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.Value(i).*[78.09/20.95;0.93/20.95];
app = Define_I(app);
%% COMPUTE PROPERTIES
app = Define_FOI(app, i);
%% PROBLEM TYPE
% Specify the problem type and the thermodynamic properties of the product
% mixture (P) required for the computations
switch app.PD.ProblemType
    case {'TP','TV'}
        TP = app.PD.TP_vector.Value;
        app.PS.strP{i} = SolveProblemTP_TV(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,TP,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
    case {'HP'}
        if i==app.C.l_phi
            % app.PS.strP{i} = SolveProblemHP_test(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value); % Newton-Raphson-worse convergence!!
            app.PS.strP{i} = SolveProblemHP(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        else
            app.PS.strP{i} = SolveProblemHP_EV_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        end
    case {'EV'}
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemEV(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        else
            app.PS.strP{i} = SolveProblemHP_EV_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        end
    case 'SP'
        pP = app.PD.pP_vector.Value(i);
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemSP(app.PS.strR{i},app.PD.phi.Value(i),pP,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        else
            app.PS.strP{i} = SolveProblemSP_fast(app.PS.strR{i},app.PD.phi.Value(i),pP,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        end
    case 'SV'
        vP_vR = app.PD.vP_vR_vector.Value(i);
        app.PS.strP{i} = SolveProblemSV(app.PS.strR{i},app.PD.phi.Value(i),vP_vR*app.PS.strR{i}.v,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
    case {'SHOCK_I','SHOCK_R'}
        u1 = app.PD.u1_vector.Value(i);
        if strcmp(app.PD.ProblemType,'SHOCK_I')
            if i==app.C.l_phi
                [app.PS.strR{i},app.PS.strP{i}] = SolveProblemSHOCK_I(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
            else
                [app.PS.strR{i},app.PS.strP{i}] = SolveProblemSHOCK_I_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
            end
        else
            if i==app.C.l_phi
                [app.PS.strR{i},app.PS.str2{i},app.PS.strP{i}] = SolveProblemSHOCK_R(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
            else
                [app.PS.strR{i},app.PS.str2{i},app.PS.strP{i}] = SolveProblemSHOCK_R_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.PS.str2{i+1},app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
            end
        end
    case 'DET'
%         app.PD.TR.Value = app.PD.TR_vector.Value(i);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_main(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
        [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_main_2(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_hybrid_main(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_Jouget_main(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp);
    case 'DET_OVERDRIVEN'
        [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_OVERDRIVEN(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,app.strThProp,app.PD.overdriven);
end
%% DISPLAY RESULTS
if ~strcmp(app.PD.ProblemType,'SHOCK_R') && ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
    displayresults(app.PS.strR{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.NameSpecies);
elseif ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
    displayresults(app.PS.strR{i},app.PS.str2{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.NameSpecies); % Display all results SHOCK_R
end
end
disp('TIME:')
toc
%% DISPLAY MOLAR FRACTION VS EQUIVALENCE RATIO
% app.M.display_species = {};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','N2','O2'};
app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb'};
% app.M.display_species = {'CN','H','Cbgrb','C2H2_acetylene','HCN','CO','C2N2','HNC','H2','N2','H2O','CO2','O2'};
% app.M.display_species = {'Cbgrb','CO2','CO','HCN','H2','OH','H2O','O2','N2'};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','N2','O2','OH','H','O','HO2','NO','HCO','CH4','CH3','HO2',...
%     'NO2','NH3','NH2','N','HCN','CH'};
closing(app,app.PS.strP,app.PD.phi.Value,app.M.display_species,app.Misc.timer_0,app.S.NameSpecies,app.C.mintol_display,app.PD.ProblemType);
%% EXCEL I/O

% delete(filename);
if app.Misc.save_Excel
    ExcelOutputMatrix_R = FormattedOutput_test([],app.PD.phi.Value,app.PS.strR,app.S.NameSpecies);
    ExcelOutputMatrix_P = FormattedOutput_test([],app.PD.phi.Value,app.PS.strP,app.S.NameSpecies);
    [STATUS1,MESSAGE1]  = xlswrite(strcat(app.C.filename,'_R.xls'),ExcelOutputMatrix_R);
    [STATUS2,MESSAGE2]  = xlswrite(strcat(app.C.filename,'_P.xls'),ExcelOutputMatrix_P);
end
% [NUM,TXT,RAW] = xlsread(filename);