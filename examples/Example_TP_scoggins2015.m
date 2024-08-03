% -------------------------------------------------------------------------
% EXAMPLE: TP_scoggins2015
%
% Compute equilibrium composition for a defined set of temperatures
% (200, 5000) [K] at atmospheric pressure (1.01325 bar) for a Si-C6H5OH_phenol
% mixture at standard conditions, and considering a set of 111 species 
% (7 in condensed phase)
% 
% This example is obtained from [1]
% 
% [1] Scoggins, J. B., & Magin, T. E. (2015). Gibbs function continuation
%     for linearly constrained multiphase equilibria. Combustion and Flame,
%     162(12), 4514-4522.
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update April 02 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

% Definitions
listSpecies = {'C','C2','C2H','C10H8_naphthale','C2H2_acetylene',...
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
    'H2O','H2ObLb','H2Obcrb','H2O2','HCCO','HCHO_formaldehy','HCN',...
    'HCO','HCOOH','HNC','HNCO','HNO','HO2','N','N2','N2O','NCO',...
    'NH','NH2','NH2OH','NH3','NO','NO2','O','O2','OCCN','OH','Si',...
    'Si2','Si2C','Si3','SiC','SiC2','SiH','SiH2','SiH3',...
    'SiH4','SiO','SiO2','SiCbbb','SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb'};

displaySpecies = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2',...
    'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb', 'H2O',...
    'SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb','H2ObLb','H2Obcrb'};

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, listSpecies);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'Si', 'C6H5OH_phenol'}, [1, 9]);

% Define properties
mixArray = setProperties(mix, 'temperature', 200:10:5000, 'pressure', 1 * 1.01325);

% Initialize solver
solver = EquilibriumSolver('problemType', 'TP');

% Solve problem
solver.solveArray(mixArray);

% Plot molar fractions
plotComposition(mixArray(1), mixArray, 'T', 'Xi', 'displaySpecies', displaySpecies, 'mintol', 1e-3);