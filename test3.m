% Example
% Import the necessary classes
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*

% Definitions
temperature = 200:10:5000;
N = length(temperature);
listSpecies = {'SiO'	'SiO2'	'SiO2ba_qzb'	'SiO2bb_qzb'	'SiO2bb_crtb'	'SiO2bLb'	'SiC'	'SiC2'	'Si2C'	'SiCbbb'	'SiCbLb'	'SiH'	'SiH2'	'SiH3'	'SiH4'	'Si'	'Si2'	'Si3'	'Sibcrb'	'SibLb'	'CH2OH'	'CH3O'	'CH3OH'	'CH3OOH'	'COOH'	'CH2CO_ketene'	'ObCHb2O'	'HObCOb2OH'	'CH3CO_acetyl'	'C2H4O_ethylen_o'	'CH3CHO_ethanal'	'CH3COOH'	'OHCH2COOH'	'C2H5OH'	'CH3OCH3'	'CH3O2CH3'	'C3H6O_propylox'	'C3H6O_acetone'	'C3H6O_propanal'	'C3H8O_1propanol'	'C3H8O_2propanol'	'bCH3COOHb2'	'C6H5O_phenoxy'	'C6H5OH_phenol'	'C7H8O_cresol_mx'	'HCO'	'HCCO'	'HCHO_formaldehy'	'HCOOH'	'bHCOOHb2'	'CH3OHbLb'	'C2H4ObLb_ethyle'	'C2H5OHbLb'	'CO'	'CO2'	'C2O'	'C3O2'	'HO2'	'H2O'	'H2O2'	'OH'	'H2Obcrb'	'H2ObLb'	'H2O2bLb'	'O'	'O2'	'O3'	'O2bLb'	'O3bLb'	'CH'	'CH2'	'CH3'	'CH4'	'C2H'	'C2H2_acetylene'	'C2H2_vinylidene'	'C2H3_vinyl'	'C2H4'	'C2H5'	'C2H6'	'C3H3_1_propynl'	'C3H3_2_propynl'	'C3H4_allene'	'C3H4_propyne'	'C3H4_cyclominus'	'C3H5_allyl'	'C3H6_propylene'	'C3H6_cyclominus'	'C3H7_n_propyl'	'C3H7_i_propyl'	'C3H8'	'C4H2_butadiyne'	'C4H4_1_3_cyclominus'	'C4H6_butadiene'	'C4H6_1butyne'	'C4H6_2butyne'	'C4H6_cyclominus'	'C4H8_1_butene'	'C4H8_cis2_buten'	'C4H8_tr2_butene'	'C4H8_isobutene'	'C4H8_cyclominus'	'C4H9_n_butyl'	'C4H9_i_butyl'	'C4H9_s_butyl'	'C4H9_t_butyl'	'C4H10_n_butane'	'C4H10_isobutane'	'C5H6_1_3cyclominus'	'C5H8_cyclominus'	'C5H10_1_pentene'	'C5H10_cyclominus'	'C5H11_pentyl'	'C5H11_t_pentyl'	'C5H12_n_pentane'	'C5H12_i_pentane'	'CH3CbCH3b2CH3'	'C6H2'	'C6H5_phenyl'	'C6H6'	'C6H10_cyclominus'	'C6H12_1_hexene'	'C6H12_cyclominus'	'C6H13_n_hexyl'	'C6H14_n_hexane'	'C7H7_benzyl'	'C7H8'	'C7H14_1_heptene'	'C7H15_n_heptyl'	'C7H16_n_heptane'	'C7H16_2_methylh'	'C8H8_styrene'	'C8H10_ethylbenz'	'C8H16_1_octene'	'C8H17_n_octyl'	'C8H18_n_octane'	'C8H18_isooctane'	'C9H19_n_nonyl'	'C10H8_naphthale'	'C10H21_n_decyl'	'C12H9_o_bipheny'	'C12H10_biphenyl'	'bCH2bxbcrb'	'CH4bLb'	'C2H2bLb_acetyle'	'C2H4bLb'	'C2H6bLb'	'C3H6bLb_propyle'	'C3H8bLb'	'C4H8bLb_1_buten'	'C4H10bLb_n_buta'	'C4H10bLb_isobut'	'C5H12bLb_n_pent'	'C6H6bLb'	'C6H14bLb_n_hexa'	'C7H8bLb'	'C7H16bLb_n_hept'	'C8H18bLb_n_octa'	'C8H18bLb_isooct'	'JP_4'	'JP_5'	'JP_10bLb'	'JP_10bgb'	'Jet_AbLb'	'Jet_Abgb'	'RP_1'	'C2H3VinylRadi'	'C7H7C6H5CH2B'	'C9H8Benzene_1minus'	'C'	'C2'	'C3'	'C4'	'C5'	'Cbgrb'	'H'	'H2'	'H2bLb'};
prefixDataName = 'C6H5OH_phenol_and_Si';
filename = {strcat(prefixDataName, '_TP1.out')};
display_species = {'H2', 'H', 'CH4', 'C2H2_acetylene', 'SiC2',...
        'Si', 'C', 'C2H', 'C3', 'CO', 'CO2', 'Cbgrb', 'SiCbbb',...
        'SiO2ba_qzb','SiO2bb_qzb','SiO2bb_crtb','H2O','H2ObLb','H2Obcrb'};
mintolDisplay = 1e-3;

tic
% Load NASA's database
DB = NasaDatabase();

% Define the chemical system, i.e., the set of species that may be present in the chemical transformation at equilibrium
system = ChemicalSystem(DB, listSpecies);

% Define the initial mixture composition
mix = Mixture(system);
set(mix, {'Si', 'C6H5OH_phenol'}, [1, 9]);

% Define the thermodynamic state
mixArray = mix.setProperties('temperature', temperature, 'pressure', 1.01325);

% Create the equilibrium solver
solver = EquilibriumSolver('flag_results', false);

% Solve chemical transformation at equilibrium
solver.solveArray(mixArray);
toc

% Load results CEA 
results_CEA = data_CEA(filename, display_species);

% Validate molar fractions
[~, fig1] = plot_molar_fractions(mixArray(1), mixArray, 'T', 'Xi', 'validation', results_CEA, 'display_species', display_species, 'mintol', mintolDisplay); % , 





