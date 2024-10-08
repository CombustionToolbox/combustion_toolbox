% -------------------------------------------------------------------------
% EXAMPLE: HP PROPELLANTS RP1
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (0.5, 5) [-]
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Feb 19 2022
% -------------------------------------------------------------------------

%% INITIALIZE
LS = {'CO2','CO','H2O','H2','O2','He','Ar','Cbgrb',...
'C2H','C2H2_acetylene','C2H3_vinyl','C2H4','C2H5OH','C2H6',...
'C3','C3H3_2_propynl','C3H4_allene','C3H4_propyne',...
'C3H5_allyl','C3H6O_acetone','C3H6_propylene',...
'C3H8','C6H2','CH','CH2CO_ketene','CH2OH','CH3CHO_ethanal',...
'CH3COOH','CH4','COOH','H','HCO','HCOOH','O','RP_1','H2bLb','O2bLb',...
'COO_M','CD_M','CDH3_M','CDO_M','CD2_M','CD2O_M','CD3_M','CD4_M',...
'CD3OD_M','CH_M','CHO_M','COH_M','HOCO_M','HCOO_M','HCbOOb_M','CHO3_M',...
'CH2_M','CH2D2_M','HCHO_M','CHOH_M','CH2O_M','CH2O2_M','CH2OO_M',...
'CH3_M','CH3OD_M','CH3O_M','CH2OH_M','CH3O2_M','CH3OH_M','CH2OH2_M',...
'CH4O2_M','COO_M','CT_M','CT3_M','CT4_M','C2_M','C2D_M',...
'C2D2_M','C2D2O_M','C2D4_M','C2D4O_M','C2D6_M','C2D6O_M','C2H_M',...
'C2DH_M','HCCO_M','C2HO3_M','C2H2_M','C2H2O_M','HCC_OH_M',...
'C2H2O2_M','O_COH_CH_O_M','C2H2O4_M','CH3C_M','CH3Cquartet_M',...
'CH3CD3_M','CH3CO_M','CH2_CHO_M','C2H3O_M','C2H3O2_M','CH3CH_M',...
'C2H4O_M','C2H4O2_M','C2H4O3_M','C2H5_M','C2H5O_M','C2H5O2_M',...
'C2H4OOH_M','C2H6O_M','C2H6O2_M','COC_M','C2O_M','C3D4_M','C3D6_M',...
'C3D6O_M','C3H_M','C3H2_M','C3H3_1_propynl_M','C3H3_M','C3H3O_M',...
'C3H4_M','C3H4O_M','C3H4O2_M','C3H5_M','C3H5O_M','C3H5O2_M',...
'C3H5OH_M','C3H6O_M','C3H6O2_M','C3H6O3_M','C3H7_n_propyl_M',...
'C3H7_i_propyl_M','C3H7O_M','C3H7OO_M','C3H7O2_M','C3H8O_M',...
'C3H8O2_M','C3H8O3_M','C3O2_M','C4_M','C4D6_M','C4H_M','C4H2_M',...
'C4H3_M','C4H4_M','C4H4O_M','num_1_4_C4H4O_M','C4H4O2_M','C4H4O4_M',...
'C4H5_M','C4H5O_M','C4H5O2_M','C4H6_M','C4H6O_M','C4H6O2_M',...
'C4H6O4_M','C4H7_M','C4H7O_M','C4H7O2_M','num_1_C4H8_M','C4H8_M',...
'C4H8O_M','C4H8O2_M','C4H7OOH_M','C4H8O4_M','C4H9_M','t_C4H9_M',...
'C4H9O_M','C4H9O2_M','C4H10O_M','num_2_C4H10O_M','T_C4H10O_M',...
'n_C4H10O2_M','C4H10O2_M','C4H10O3_M','C4H10O4_M','C5_M','C5H_M',...
'C5H2_M','C5H3_M','C5H4_M','C5H4O_M','C5H4O2_M','C5H5_M','C5H4OH_M',...
'C5H5O_M','C5H5O2_M','C5H6_M','C5H5OH_M','C5H6O_M','C5H6O2_M',...
'C5H7_M','C5H7O_M','C5H8_M','C5H8O_M','C5H8O2_M','C5H8O4_M',...
'C5H9_M','C5H9O2_M','C5H10_M','C5H9OH_M','C5H10O_M','C5H10O2_M',...
'C5H10O3_M','n_C5H11_M','s_C5H11_M','t_C5H11_M','C5H11_M',...
'C5H11OH_M','C5H12O_M','C5H12O2_M','C6_M','C6H_M','C6H3_M',...
'o_C6H3_M','num_1_2_C6H4_M','num_1_3_C6H4_M','num_1_4_C6H4_M',...
'num_1_5_C6H4_M','num_1_2_3_4_5_C6H4_M','C6H4_M','C6H4O2_M',...
'C6H5_M','C6H5O_M','C6H5OO_M','C6H6_M','C6H6O_M','C6H6O2_M',...
'C6H5OOH_M','C6bOHb6_M','C6H7_M','C6H8_M','C6H8O_M','C6H8O3_M',...
'C6H8O7_M','C6H9_M','C6H10_M','C6H10O5_M','C6H11_M','C6H11O2_M',...
'C6H12_M','C6H12O_M','C6H12O2_M','C6H12O6_M','C6H12O7_M','C6H13_M',...
'C6H14_M','C6H14O_M','C6H14O2_M','C6H14O3_M','C6H14O6_M','C6T6_M',...
'C7_M','C7H_M','C7H4_M','C7H5O_M','C7H6O_M','C7H6O2_M','C7H7_M',...
'C7H7O_M','C7H7O2_M','C7H8_M','C7H8O_M','C7H8O2_M','C7H9_M','C7H10_M',...
'C7H11_M','C7H12_M','C7H13_M','C7H14_M','C7H14O_M','C7H14O2_M',...
'C7H15_M','C7H16_M','C7H16O_M','C8_M','C8H_M','C8H2_M','C8H5_M',...
'C8H6_M','C8H6O_M','C8H6O2_M','C8H7_M','C8H8_M','C8H8O_M','C8H8O2_M',...
'C8H9_M','C8H10_M','C8H12_M','C8H14_M','C8H16_M','C8H16O2_M',...
'C8H17_M','C8H18_M','C8H18O_M','C9_M','C9H_M','C9H4_M','C9H7_M',...
'C9H8_M','C9H9_M','C9H10_M','C9H10O2_M','C9H11_M','C9H12_M',...
'C9H16_M','C9H18O2_M','C9H18O6_M','C9H19_M','C9H20_M','C9H20O_M',...
'C10_M','C10D8_M','C10H_M','C10H2_M','C10H6_M','C10H7_M',...
'C10H8_M','C10H8O_M','C10H9_M','C10H10_M','C10H12O3_M','C10H14_M',...
'C10H15_M','C10H16_exo_M','C10H16_M','C10H18_M','C10H20_M',...
'C10H20O2_M','C10H21_1_M','C10H21_M','C10H22_M','C10H22O_M',...
'C10H22O4_M','C11_M','C11H_M','C11H8_M','C11H9_M','C11H10_M',...
'C11H11_M','C11H12_M','C11H20O2_M','C11H22O2_M','C11H24_M','C12_M',...
'O_C12D9_M','C12D10_M','C12H_M','C12H2_M','C12H7_M','C12H8_M',...
'C12H8O_M','C12H8O2_M','C12H9_M','C12H10_M','C12H12_M','C12H18_M',...
'C12H24_M','C12H24O2_M','C12H25_M','C12H26_M','C13H9_M','C13H10_M',...
'C13H10O_M','C13H12_M','p_C13H12_M','C13H14_M','C13H26O2_M',...
'C14H8_M','C14H10_M','C14H12_M','C14H14_M','C14H28_M','C14H28O2_M',...
'n_C14H30_M','C15H12_M','C15H14_M','C15H16O2_M','C15H30O2_M',...
'n_C15H32_M','C16H10_M','C16H29O2_M','C16H30O2_M','C16H31O2_M',...
'C16H32O2_M','C16H33_M','C16H34_M','C16H34O_M','C17H12_M',...
'C17H31O2_M','C17H32O2_M','C17H33O2_M','C17H34O2_M','C17H36_M',...
'C18H10_M','C18H12_M','C18H29O2_M','C18H30O2_M','C18H31O2_M',...
'C18H32O2_M','C18H33O2_M','C18H34_M','C18H34O2_M','C18H34O3_M',...
'C18H35O2_M','C18H36_M','C18H36O2_M','C18H36O4_M','C18H38_M',...
'C19H40_M','C20H10_M','C20H12_M','C20H32O2_M','C20H38O2_M',...
'C20H39O2_M','C20H40O2_M','C20H42_M','C21H12_M','C21H44_M',...
'C22H14_M','C22H44O2_M','C24H12_M','C24H18_M','C24H46O2_M',...
'C24H48O2_M','C25H20_M','C25H52_M','C30H10_M','C30H62_M','C32H13_M',...
'C32H14_M','C60bcrb_M','C60_M','C70bcrb_M','C70_M','OH_M','HOT_M',...
'HO2_M','HO3_M','HT_M','H2O2_M','H2O3_M','HOOOH_M'};


% LS = {'C10H8_M', 'C', 'C2_M', 'H2', 'H2O', 'C6H6_M', 'CO', 'C7H8_M', 'CH4', 'O2',...
%     'CO2', 'OH_M', 'RP_1', 'O2bLb'};

% 'HC/O2/N2 PROPELLANTS'
self = App(LS);
self.TN.tolN = 1e-20;
self.TN.tol0 = 1e-6;
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 30, 'phi', 6.811333:0.5:11.352221); % RP-1
% self = set_prop(self, 'TR', 300, 'pR', 172.35, 'phi', 6:0.5:10); % 'C11H20O2_M', 'C10H22_M', 'C12H26_M', 'C13H26O2_M'
self.PD.S_Fuel     = {'RP_1'};
self.PD.N_Fuel     = [1];
% self.PD.S_Fuel     = {'C11H20O2_M', 'C10H22_M', 'C12H26_M', 'C13H26O2_M'};
% self.PD.N_Fuel     = [0.354, 0.15, 0.183, 0.313];
self.PD.S_Oxidizer = {'O2bLb'};
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value); 
load huzel OF_T_real T_real OF_gamma_real gamma_real
self.Tcurve = griddedInterpolant(OF_T_real, T_real, 'pchip', 'linear');
self.Tcurve = griddedInterpolant([0.35, 0.5], [850, 1200], 'pchip', 'linear');
self.gammacurve = griddedInterpolant(OF_gamma_real, gamma_real, 'pchip', 'linear');
% self.eta_c = 0.9;
%% SOLVE PROBLEM
self = solve_problem(self, 'HP');
%% DISPLAY RESULTS (PLOTS)
self.C.mintol_display = 1e-6;
post_results(self);

nfrec = 3;
OF_CT = cell2vector(self.PS.strR, 'OF');
cp_CT = cell2vector(self.PS.strP, 'cP');
m_CT = cell2vector(self.PS.strP, 'mi');
cp_mass_CT = cp_CT ./ m_CT;

% Temperature
self.Misc.config.labelx = 'Mixing ratio $O/F$';
self.Misc.config.labely = 'Temperature $T$ [K]';
self.Misc.config.title = 'Turbine inlet temperature [K]';
ax = plot_figure(OF_CT, self.PS.strP, 'OF', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
plot(ax, OF_T_real(1:nfrec:end), T_real(1:nfrec:end), 'ko', 'MarkerFaceColor', 'auto');
legend_name = {'CT', 'Huzel'};
set_legends(ax, legend_name, self.Misc.config)
% Adiabatic index
self.Misc.config.labely = 'Adiabatic index $\gamma$';
self.Misc.config.title = 'Adiabatic index $\gamma$';
ax = plot_figure(OF_CT, self.PS.strP, 'OF', 'gamma', self.Misc.config, self.PD.CompleteOrIncomplete);
plot(ax, OF_gamma_real(1:nfrec:end), gamma_real(1:nfrec:end), 'ko', 'MarkerFaceColor', 'auto');
legend_name = {'CT', 'Huzel'};
set_legends(ax, legend_name, self.Misc.config)
% Heat capacity at constant pressure
self.Misc.config.labely = '$c_p$ [J/kg-K]';
self.Misc.config.title = 'Heat capacity at constant pressure [J/kg-K]';
plot_figure(OF_CT, cp_mass_CT, 'OF', 'cP');
% Mean molecular weight
plot_figure(OF_CT, self.PS.strP, 'OF', 'W');
% Efficiency in the chamber 
plot_figure(OF_CT, self.PS.strP, 'OF', 'eta_c');