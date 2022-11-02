% -------------------------------------------------------------------------
% EXAMPLE: TP
%
% Compute equilibrium composition at defined temperature (e.g., 3000 K) and
% pressure (e.g., 1.01325 bar) for lean to rich CH4-air mixtures at
% standard conditions, a set of 26 species considered and a set of
% equivalence ratios phi contained in (0.5, 5) [-]
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','He','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update July 22 2022
% -------------------------------------------------------------------------
LS = {'Ar','C','C10H8_naphthale','C2','C2H','C2H2_acetylene',...
    'C2H2_vinylidene','C2H3_vinyl','C2H4','C2H4O_ethylen_o','C2H5',...
    'C2H5OH','C2H6','C2N2','C2O','C3','C3H3_1_propynl',...
    'C3H3_2_propynl','C3H4_allene','C3H4_propyne','C3H5_allyl',...
    'C3H6O_acetone','C3H6_propylene','C3H8','C3O2','C4',...
    'C4H2_butadiyne','C4H6_1butyne','C4H6_2butyne','C4H6_butadiene',...
    'C4H8_1_butene','C4H8_cis2_buten','C4H8_isobutene',...
    'C4H8_tr2_butene','C4H9_i_butyl','C4H9_n_butyl','C4H9_s_butyl',...
    'C4H9_t_butyl','C4N2','C5','C6H2','C6H5OH_phenol','C6H5O_phenoxy',...
    'C6H5_phenyl','C6H6','C7H7_benzyl','C7H8','C8H18_isooctane',...
    'C8H8_styrene','CH','CH2','CH2CO_ketene','CH2OH','CH3',...
    'CH3CHO_ethanal','CH3CN','CH3COOH','CH3CO_acetyl','CH3O',...
    'CH3OCH3','CH3OH','CH4','CN','CO','CO2','COOH','Cbgrb','H','H2',...
    'H2O','H2O2','H2ObLb','H2Obcrb','HCCO','HCHO_formaldehy','HCN',...
    'HCO','HCOOH','HNC','HNCO','HNO','HO2','He','N','N2','N2O','NCO',...
    'NH','NH2','NH2OH','NH3','NO','NO2','O','O2','OCCN','OH','Si',...
    'Si2','Si2C','Si3','SiC','SiC2','SiCbbb','SiH','SiH2','SiH3',...
    'SiH4','SiO','SiO2','SiO2ba_qzb','SiO2bb_crtb','SiO2bb_qzb'};
%% INITIALIZE
self = App(LS);
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1.01325);
self.PD.S_Fuel = {'C6H5OH_phenol', 'Si'};
self.PD.N_Fuel = [9, 1];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 2000:50:5000);
self.Misc.display_species = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2', 'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb'};
self.C.mintol_display = 1e-4;
self.TN.tolN = 1e-16;
%% SOLVE PROBLEM
self = solve_problem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
post_results(self);