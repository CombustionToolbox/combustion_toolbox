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
%
% Authors:
%
% * Alberto Cuadra Lara, Universidad Carlos III de Madrid (UC3M)
% * Marcos Vera Coello,  Universidad Carlos III de Madrid (UC3M)
%
% Last update: 15-Jan-2020 10:57
%% LOAD DATABASES AND GLOBAL CONSTANTS
addpath(genpath(pwd));
[app,strThProp,strMaster] = Initialize(); 
app.Misc.save_Excel = false;
%% REACTION: COMPLETE OR INCOMPLETE
% Specify the type of the reaction: complete or incomplete (dissociation)

% app.PD.CompleteOrIncomplete = 'complete';
app.PD.CompleteOrIncomplete = 'incomplete'; 
app.TN.factor_c = 0.8; % factor_c = 1 (default).
%% MINORS PRODUCTS
% Specify the minority products to be considered in the product mixture (P) in
% addition to the major species (CO2, CO, H2O, H2, O2, N2, C(gr))
% Also, He and Ar are always included in the database.

% app.M.minor_products = {'O','OH','H','NO','N'};
% app.M.minor_products = {'O2', 'N2', 'CH4', 'N', 'CO2', 'H2O', 'CO', 'H2', 'O',...
%     'H', 'OH', 'NO', 'NO2', 'N2O', 'He', 'Ar', 'C3H8', 'H2O2', 'HCN', 'HNC',...
%     'CH3','CH','C2H2_acetylene','C6H6'};
%%
% *airNASA9ions.cti + others*

% app.M.minor_products = {'O2','N2','O','N','NO','C','C2','CO','CO2','CN',...
%     'Ar','CH4','H2O','H2','H','OH','He'};
%%
% *gri30-x.cti except 'CH2(s)' + others*

% app.M.minor_products = {'H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2',...
%     'CH', 'CH2', 'CH3', 'CH4', 'CO', 'CO2', 'HCO',...
%     'CH2OH', 'CH3O', 'CH3OH', 'C2H', 'C2H4',...
%     'C2H5', 'C2H6', 'HCCO', 'N', 'NH', 'NH2', 'NH3',...
%     'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN',...
%     'NCO', 'N2', 'Ar', 'C3H8','C2','C2H2_acetylene','C6H6',...
%     'C8H18_isooctane','C2H5OH','He','HNC'};
%%
% HC/O2/N2 EXTENDED

% app.M.minor_products = {'OH','H','O','HO2','NO','HCO','CH4','CH3','HO2',...
%     'NO2','NH3','NH2','N','HCN','CN','N2O','C2','CH'};
%%
% load strThProp_HC_47.mat
% SOOT FORMATION 
app.M.minor_products = {'H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2',...
    'CH', 'CH2', 'CH3', 'CH4', 'CO', 'CO2', 'HCO',...
    'CH2OH', 'CH3O', 'CH3OH', 'C2H', 'C2H4',...
    'C2H5', 'C2H6', 'HCCO', 'N', 'NH', 'NH2', 'NH3',...
    'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN',...
    'NCO', 'N2', 'Ar', 'C3H8','C2','C2H2_acetylene','C6H6',...
    'C8H18_isooctane','C2H5OH','He','HNC','HNCO','NH2OH',...
    };
% app.M.minor_products = fieldnames(strThProp)';
% SOOT FORMATION WITHOUT CH4 CONSIDERED AS MINOR SPECIE
% app.M.minor_products = {'H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2',...
%     'CH', 'CH2', 'CH3', 'CO', 'CO2', 'HCO',...
%     'CH2OH', 'CH3O', 'CH3OH', 'C2H', 'C2H4',...
%     'C2H5', 'C2H6', 'HCCO', 'N', 'NH', 'NH2', 'NH3',...
%     'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN',...
%     'NCO', 'N2', 'Ar', 'C3H8','C2','C2H2_acetylene','C6H6',...
%     'C8H18_isooctane','C2H5OH','He','HNC','HNCO','NH2OH'};
%%
% *NASA ALL: CH4 + 2O2 + 7.52N2*

% app.M.minor_products = {'CH3','CH4','CN','CO2','C2H','CH2CO_ketene',...
%     'C2H3_vinyl','C2H4','CH3COOH','C2H6','CH3OCH3','CNC','C2O',...
%     'C3H3_2_propynl','C3H4_cyclominus','C3H6_cyclominus','C3H6O_propanal',...
%     'C3H8','CNCOCN','C4H2_butadiyne','C4H6_1butyne','C4H8_1_butene',...
%     'C4H8_isobutene','C4H9_n_butyl','C4H9_t_butyl','C4N2','C5H8_cyclominus',...
%     'C5H11_pentyl','C5H12_i_pentane','C6H5_phenyl','C6H5OH_phenol',...
%     'C6H12_cyclominus','C7H7_benzyl','C7H14_1_heptene','C7H16_2_methylh',...
%     'C8H16_1_octene','C8H18_isooctane','C10H21_n_decyl','H','HCCN','HNCO',...
%     'HNO3','HCHO_formaldehy','H2O2','NCO','NH3','NO2','N2O','NH2NO2','N2O4',...
%     'N3H','O2','C2H5OHbLb','C6H6bLb','H2ObLb','CH','CH2OH','CH3OH',...
%     'CNN','COOH','C2H2_acetylene','ObCHb2O','CH3CN','C2H4O_ethylen_o',...
%     'OHCH2COOH','CH3N2CH3','CH3O2CH3','OCCN','C3','C3H4_allene','C3H5_allyl',...
%     'C3H6O_propylox','C3H7_n_propyl','C3H8O_1propanol','C3O2',...
%     'C4H4_1_3_cyclominus','C4H6_2butyne','C4H8_cis2_buten','C4H8_cyclominus',...
%     'C4H9_i_butyl','C4H10_n_butane','C5','C5H10_1_pentene','C5H11_t_pentyl',...
%     'CH3CbCH3b2CH3','C6H5O_phenoxy','C6H10_cyclominus','C6H13_n_hexyl','C7H8',...
%     'C7H15_n_heptyl','C8H8_styrene','C8H17_n_octyl','C9H19_n_nonyl','C12H9_o_bipheny',...
%     'HCN','HCCO','HNO','HO2','HCOOH','NH','NH2OH','NO3','NCN','N2H4',...
%     'N2O5','O','O3','N2H4bLb','CH2',...
%     'CH3O','CH3OOH','CO','C2','C2H2_vinylidene','HObCOb2OH','CH3CO_acetyl',...
%     'CH3CHO_ethanal','C2H5','C2H5OH','CCN','C2N2','C3H3_1_propynl','C3H4_propyne',...
%     'C3H6_propylene','C3H6O_acetone','C3H7_i_propyl','C3H8O_2propanol',...
%     'C4','C4H6_butadiene','C4H6_cyclominus','C4H8_tr2_butene',...
%     'C4H9_s_butyl','C4H10_isobutane','C5H6_1_3cyclominus','C5H10_cyclominus',...
%     'C5H12_n_pentane','C6H2','C6H6','C6H12_1_hexene','C6H14_n_hexane',...
%     'C7H8O_cresol_mx','C7H16_n_heptane','C8H10_ethylbenz','C8H18_n_octane',...
%     'C10H8_naphthale','C12H10_biphenyl','HCO','HNC','HNO2','H2','H2O',...
%     'N','NH2','NO','N2','N2H2','N2O3','N3','OH','CH3OHbLb','C6H5NH2bLb','He','Ar','C'};
%%%%%%%%%%%%%%%
% app.M.minor_products = {'CH3','CN','CO2','C2H','CH2CO_ketene',...
%     'C2H3_vinyl','C2H4','CH3COOH','C2H6','CH3OCH3','CNC','C2O',...
%     'C3H3_2_propynl','C3H4_cyclominus','C3H6_cyclominus','C3H6O_propanal',...
%     'C3H8','CNCOCN','C4H2_butadiyne','C4H6_1butyne','C4H8_1_butene',...
%     'C4H8_isobutene','C4H9_n_butyl','C4H9_t_butyl','C4N2','C5H8_cyclominus',...
%     'C5H11_pentyl','C5H12_i_pentane','C6H5_phenyl','C6H5OH_phenol',...
%     'C6H12_cyclominus','C7H7_benzyl','C7H14_1_heptene','C7H16_2_methylh',...
%     'C8H16_1_octene','C8H18_isooctane','C10H21_n_decyl','H','HCCN','HNCO',...
%     'HNO3','HCHO_formaldehy','H2O2','NCO','NH3','NO2','N2O','NH2NO2','N2O4',...
%     'N3H','O2','CH','CH2OH','CH3OH',...
%     'CNN','COOH','C2H2_acetylene','ObCHb2O','CH3CN','C2H4O_ethylen_o',...
%     'OHCH2COOH','CH3N2CH3','CH3O2CH3','OCCN','C3','C3H4_allene','C3H5_allyl',...
%     'C3H6O_propylox','C3H7_n_propyl','C3H8O_1propanol','C3O2',...
%     'C4H4_1_3_cyclominus','C4H6_2butyne','C4H8_cis2_buten','C4H8_cyclominus',...
%     'C4H9_i_butyl','C4H10_n_butane','C5','C5H10_1_pentene','C5H11_t_pentyl',...
%     'CH3CbCH3b2CH3','C6H5O_phenoxy','C6H10_cyclominus','C6H13_n_hexyl','C7H8',...
%     'C7H15_n_heptyl','C8H8_styrene','C8H17_n_octyl','C9H19_n_nonyl','C12H9_o_bipheny',...
%     'HCN','HCCO','HNO','HO2','HCOOH','NH','NH2OH','NO3','NCN','N2H4',...
%     'N2O5','O','O3','N2H4bLb','CH2',...
%     'CH3O','CH3OOH','CO','C2','C2H2_vinylidene','HObCOb2OH','CH3CO_acetyl',...
%     'CH3CHO_ethanal','C2H5','C2H5OH','CCN','C2N2','C3H3_1_propynl','C3H4_propyne',...
%     'C3H6_propylene','C3H6O_acetone','C3H7_i_propyl','C3H8O_2propanol',...
%     'C4','C4H6_butadiene','C4H6_cyclominus','C4H8_tr2_butene',...
%     'C4H9_s_butyl','C4H10_isobutane','C5H6_1_3cyclominus','C5H10_cyclominus',...
%     'C5H12_n_pentane','C6H2','C6H6','C6H12_1_hexene','C6H14_n_hexane',...
%     'C7H8O_cresol_mx','C7H16_n_heptane','C8H10_ethylbenz','C8H18_n_octane',...
%     'C10H8_naphthale','C12H10_biphenyl','HCO','HNC','HNO2','H2','H2O',...
%     'N','NH2','NO','N2','N2H2','N2O3','N3','OH','C6H5NH2bLb','He','Ar'};
%%%%%%%%%%%%%%%%%
% ,'H2ObLb' ,'C2H5OHbLb' ,'CH3OHbLb','C6H6bLb'
%%
%  *NASA : CH4 + 2O2 + 7.52N2*

% app.M.minor_products = {'C','CN','CO2','H','O2','CH','C3','C5','NH','O','CO',...
%     'C2','C4','H2','N','NO','N2','OH','CH4','H2O','He','Ar'};
%%
% *SHOCK NASA O2+N2 + OTHERS*

% app.M.minor_products = {'O2','N2','O','O3','N','NO','NO2','NO3','N2O3','N2O4','N3','C','CO','CO2',...
%     'Ar','CH4','CH3','CH','H2O','H2','H','He'};
%%
% *HYDROGEN*

%  app.M.minor_products = {'H','HNO3','H2O','NH','NH2OH','NO3','N2H2','N2O3','N3','OH','HNO2',...
%     'H2','N','NH3','NO2','N2O','N2H4','N2O5','O','O3','He','Ar','CO2','CO','O2','N2','HO2','NH2','H2O2'};
%% CHECK SPECIES MINOR PRODUCTS --> SELECTED DATABASE (strThProp)
% Checks if any specie of the app.M.minor_products considered are not included
% in the selected database strThProp. In case some is missing compute it.
% IF THE MINOR-SPECIE IS NOT DECLARE AS REACTANT THE THEORY WILL FAIL!!
[strThProp,app.E,app.S,app.M,app.C] = check_database(strMaster,strThProp,app.E,app.S,app.M,app.C);
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
app.PD.pR.Value = 1;
% app.PD.phi.Value = 1*ones(1,length(app.PD.TR.vector.Value));
app.PD.phi.Value = 0.1:0.1:2;
% app.PD.phi.Value = 2.5:0.001:2.64;
% app.PD.phi.Value = 0.1:0.1:0.3;
% app.PD.phi.Value = 2;
%% INITIALIZATION
[app.E,app.S,app.M,app.C,app.Misc,Problem_selected] = Initialize_2(app.E,app.S,app.M,app.C,app.Misc);
%% PROBLEM TYPE
switch Problem_selected
    case 'TP' % * TP: Equilibrium composition at defined T and p
        app.PD.ProblemType = 'TP';
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
for i=app.C.l_phi:-1:1 % Evading preallocate struct
% waitbar(1-i/app.C.l_phi,f,strcat('Case ',sprintf(' %d',app.C.l_phi-i+1),' -',sprintf(' %d',app.C.l_phi)));
% if getappdata(f,'canceling'), break, end
%% DEFINE FUEL
% Define the fuel blend by specifying the fuel(s), the number of moles and
% the temperature T [K]
% app.PD.TR.Value = app.PD.TR.vector.Value(i);
% app.PD.R_Fuel = 0; app.PD.phi_t = 1; app.PD.Fuel.x = 0; app.PD.Fuel.eps = 1e-1;
% app.PD.S_Fuel = {'CH4','C2H6','C3H8'}; app.PD.N_Fuel = [0.85;0.1;0.05]; 
app.PD.S_Fuel = {'CH4'}; app.PD.N_Fuel = 1;
% app.PD.S_Fuel = {'H2'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C2H2_acetylene'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C3H8'}; app.PD.N_Fuel = 1;    
% app.PD.S_Fuel = {'C6H6'}; app.PD.N_Fuel = 1; 
% app.PD.S_Fuel = {'C2H6'}; app.PD.N_Fuel = 1;
% app.PD.S_Fuel = {'C2H5OH'}; app.PD.N_Fuel = 1;

if ~isfield(app.PD,'R_Fuel')
    app.PD.R_Fuel = SetSpecies(app.C.M0.Value,app.PD.S_Fuel,app.PD.N_Fuel,app.PD.TR.Value,find_idx(app.PD.S_Fuel,app.S.NameSpecies),strThProp);
    app.PS.strR_Fuel{i} = ComputeProperties(app.C.A0.Value,app.PD.R_Fuel,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H); app.PS.strR{i}.phi = app.PD.phi.Value(i);
    app.PD.Fuel.x = app.PS.strR_Fuel{i}.NatomE(app.E.ind_C);
    app.PD.Fuel.y = app.PS.strR_Fuel{i}.NatomE(app.E.ind_H);
    app.PD.Fuel.z = app.PS.strR_Fuel{i}.NatomE(app.E.ind_O);
    app.PD.Fuel.eps = 0;
    app.PD.phi_t = app.PD.Fuel.x+app.PD.Fuel.y/4-app.PD.Fuel.z/2;
end
% DG0 = (species_g0_new('CO',TP,strThProp)-species_g0_new('CO2',TP,strThProp))*1000;
% k_c = exp(-DG0/(R0*TP));
% eps = 1/k_c;
% eps = 1e-3;
% app.PD.phi_c = -(((-1+eps)*(4*x+y-2*z))/(2*(1+eps)*x+eps^2*y-2*z+2*eps*z)); % C_x H_y O_z
% app.PD.phi_c = 2/(x)*(x+y/4);       % C_x H_y
%% DEFINE OXIDIZER
% Specify the oxidizer, either i) pure O2, specifying the number of moles
% of 'O2' or the equivalence ratio 'O2Phi', or ii) air, specifying the
% equivalence ratio 'AirPhi'. When specifying the equivalence ratio one
% must provide R_fuel as an input so that the required ammount of oxygen
% for stoichiometric fuel-O2 combustion can be evaluated first

% app.PD.R_Oxidizer = 0;
% app.PD.R_Oxidizer = SetSpecies('O2',7.5/app.PD.phi.Value(i),'T',app.PD.TR.Value);
% app.PD.R_Oxidizer = SetSpecies('Air',1/app.PD.phi.Value(i),'T',app.PD.TR.Value);
% app.PD.R_Oxidizer = SetSpecies('O2Phi',app.PD.R_Fuel,app.PD.phi.Value(i),'T',app.PD.TR.Value);
% app.PD.R_Oxidizer = SetSpecies('AirPhi',app.PD.R_Fuel,app.PD.phi.Value(i),'T',app.PD.TR.Value);
app.PD.S_Oxidizer = {'O2'}; app.PD.N_Oxidizer = app.PD.phi_t/app.PD.phi.Value(i);
app.PD.R_Oxidizer = SetSpecies(app.C.M0.Value,app.PD.S_Oxidizer,app.PD.N_Oxidizer,app.PD.TR.Value,find_idx(app.PD.S_Oxidizer,app.S.NameSpecies),strThProp);
%% DEFINE DILUENTS/INERTS
% Specify diluents/inert species, such as CO2, H2O, N2, He, Ar, etc.

% app.PD.R_Inert = 0; % case without inert gases
app.PD.S_Inert = {'N2'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.Value(i)*79/21;
% app.PD.S_Inert = {'N2','Ar'}; app.PD.N_Inert = app.PD.phi_t/app.PD.phi.Value(i).*[78.09/20.95;0.93/20.95];
app.PD.R_Inert = SetSpecies(app.C.M0.Value,app.PD.S_Inert,app.PD.N_Inert,app.PD.TR.Value,find_idx(app.PD.S_Inert,app.S.NameSpecies),strThProp);
%% COMPUTE PROPERTIES
% COMPUTE PROPERTIES OF THE REACTIVES FOR THE GIVEN CONDITIONS
R = app.PD.R_Fuel+app.PD.R_Oxidizer+app.PD.R_Inert;
app.PS.strR{i} = ComputeProperties(app.C.A0.Value,R,app.PD.pR.Value,app.PD.TR.Value,app.E.ind_C,app.E.ind_H); app.PS.strR{i}.phi = app.PD.phi.Value(i);

% FO = sum(app.PD.R_Fuel(:,1))/sum(app.PD.R_Oxidizer(:,1));
% FO_st = sum(app.PD.R_Fuel(:,1))/(app.PD.Fuel.x+app.PD.Fuel.y/4-app.PD.Fuel.z/2);
%% PROBLEM TYPE
% Specify the problem type and the thermodynamic properties of the product
% mixture (P) required for the computations
switch app.PD.ProblemType
    case {'TP','TV'}
        TP = app.PD.TP_vector.Value;
        app.PS.strP{i} = SolveProblemTP_TV(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,TP,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
    case {'HP'}
        if i==app.C.l_phi
            % app.PS.strP{i} = SolveProblemHP_test(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value); % Newton-Raphson-worse convergence!!
            app.PS.strP{i} = SolveProblemHP(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        else
            app.PS.strP{i} = SolveProblemHP_EV_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        end
    case {'EV'}
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemEV(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        else
            app.PS.strP{i} = SolveProblemHP_EV_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        end
    case 'SP'
        pP = app.PD.pP_vector.Value(i);
        if i==app.C.l_phi
            app.PS.strP{i} = SolveProblemSP(app.PS.strR{i},app.PD.phi.Value(i),pP,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        else
            app.PS.strP{i} = SolveProblemSP_fast(app.PS.strR{i},app.PD.phi.Value(i),pP,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        end
    case 'SV'
        vP_vR = app.PD.vP_vR_vector.Value(i);
        app.PS.strP{i} = SolveProblemSV(app.PS.strR{i},app.PD.phi.Value(i),vP_vR*app.PS.strR{i}.v,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
    case {'SHOCK_I','SHOCK_R'}
        u1 = app.PD.u1_vector.Value(i);
        if strcmp(app.PD.ProblemType,'SHOCK_I')
            if i==app.C.l_phi
                [app.PS.strR{i},app.PS.strP{i}] = SolveProblemSHOCK_I(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
            else
                [app.PS.strR{i},app.PS.strP{i}] = SolveProblemSHOCK_I_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
            end
        else
            if i==app.C.l_phi
                [app.PS.strR{i},app.PS.str2{i},app.PS.strP{i}] = SolveProblemSHOCK_R(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
            else
                [app.PS.strR{i},app.PS.str2{i},app.PS.strP{i}] = SolveProblemSHOCK_R_fast(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,u1,app.PS.str2{i+1},app.PS.strP{i+1},app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
            end
        end
    case 'DET'
%         app.PD.TR.Value = app.PD.TR_vector.Value(i);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_main(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
        [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_main_2(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_hybrid_main(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
%         [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_Jouget_main(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp);
    case 'DET_OVERDRIVEN'
        [app.PS.strR{i},app.PS.strP{i},app.TN.guess] = SolveProblemDET_OVERDRIVEN(app.PS.strR{i},app.PD.phi.Value(i),app.PD.pR.Value,app.PD.TR.Value,app.TN.guess,app.E,app.S,app.C,app.M,app.PD,app.TN,strThProp,app.PD.overdriven);
end
%% DISPLAY RESULTS
if ~strcmp(app.PD.ProblemType,'SHOCK_R') && ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
    displayresults(app.PS.strR{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.NameSpecies);
elseif ~strcmp(app.PD.ProblemType,'DET_OVERDRIVEN')
    displayresults(app.PS.strR{i},app.PS.str2{i},app.PS.strP{i},app.PD.ProblemType,app.C.mintol_display,app.S.NameSpecies); % Display all results SHOCK_R
end
end
%% DISPLAY MOLAR FRACTION VS EQUIVALENCE RATIO
app.M.display_species = {};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','N2','O2'};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb'};
% app.M.display_species = {'CN','H','Cbgrb','C2H2_acetylene','HCN','CO','C2N2','HNC','H2','N2','H2O','CO2','O2'};
app.M.display_species = {'Cbgrb','CO2','CO','HCN','H2','OH','H2O','O2','N2'};
% app.M.display_species = {'CO','CO2','H','HO2','H2','H2O','N2','O2','OH','H','O','HO2','NO','HCO','CH4','CH3','HO2',...
%     'NO2','NH3','NH2','N','HCN','CH'};
closing(app.PS.strP,app.PD.phi.Value,app.M.display_species,app.Misc.timer_0,app.S.NameSpecies,app.C.mintol_display,app.PD.ProblemType);
%% EXCEL I/O

% delete(filename);
if app.Misc.save_Excel
    ExcelOutputMatrix_R = FormattedOutput_test([],app.PD.phi.Value,app.PS.strR,app.S.NameSpecies);
    ExcelOutputMatrix_P = FormattedOutput_test([],app.PD.phi.Value,app.PS.strP,app.S.NameSpecies);
    [STATUS1,MESSAGE1]  = xlswrite(strcat(app.C.filename,'_R.xls'),ExcelOutputMatrix_R);
    [STATUS2,MESSAGE2]  = xlswrite(strcat(app.C.filename,'_P.xls'),ExcelOutputMatrix_P);
end
% [NUM,TXT,RAW] = xlsread(filename);