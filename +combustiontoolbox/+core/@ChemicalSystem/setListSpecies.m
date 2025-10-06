function [obj, listSpecies, listSpeciesFormula] = setListSpecies(obj, varargin)
    % Set list of species in the mixture (products)
    %
    % Predefined list of species:
    %     * SOOT FORMATION
    %     * COMPLETE
    %     * HC/O2/N2 EXTENDED
    %     * SOOT FORMATION EXTENDED
    %     * AIR, DISSOCIATED AIR
    %     * AIR IONS, AIR_IONS
    %     * IDEAL_AIR, AIR_IDEAL
    %     * HYDROGEN
    %     * HYDROGEN_L, HYDROGEN (L)
    %     * HC/O2/N2 PROPELLANTS
    %     * SI/HC/O2/N2 PROPELLANTS
    %
    % Args:
    %     obj (ChemicalSystem): Chemical system object
    %
    % Optional Args:
    %     * listSpecies (cell): Name list species / list of species
    %     * phi (float): Equivalence ratio
    %     * phi_c (float): Equivalence ratio in which theoretically appears soot
    %
    % Returns:
    %     Tuple containing
    %
    %     * obj (ChemicalSystem): Chemical system object
    %     * listSpecies (cell): List of species
    %     * listSpeciesFormula (cell): List with chemical formula of species
    %
    % Examples:
    %     * [obj, listSpecies, listSpeciesFormula] = setListSpecies(chemicalSystem, 'soot formation');
    %     * [obj, listSpecies, listSpeciesFormula] = setListSpecies(chemicalSystem, 'soot formation');
    %     * [obj, listSpecies, listSpeciesFormula] = setListSpecies(chemicalSystem, 'complete', 1.5, 2.5);
    
    % Definitions
    database = obj.database;

    % Initialization
    if isempty(varargin)
        listSpecies = [];
        return
    end
    
    % Get list of species
    [obj.listSpecies, obj.FLAG_COMPLETE] = getListSpecies(obj, database, varargin{:});
    
    % Check that there are not repeated species
    obj.listSpecies = unique(obj.listSpecies, 'stable');

    % Set list of species
    listSpecies = obj.listSpecies;
    listSpeciesFormula = getFormula(obj.listSpecies, database);
    
    % Assign values
    obj.listSpeciesFormula = listSpeciesFormula;

    % Check if the list of species contains ions
    if any(obj.getIndexIons(obj.listSpecies))
        obj.FLAG_ION = true;
    end

end

% SUB-PASS FUNCTIONS
function [listSpecies, FLAG_COMPLETE] = getListSpecies(obj, database, varargin)
    % Get list of species
    
    % Default values
    FLAG_COMPLETE = false;

    % Unpack additional inputs
    listSpecies = varargin{1};

    % Check input
    if iscell(listSpecies)
        return
    elseif any(combustiontoolbox.utils.findIndex(database.listSpecies, listSpecies))
        listSpecies = {listSpecies};
        return
    end
    
    % Check default list of species
    switch upper(listSpecies)
        case {'COMPLETE', 'COMPLETE REACTION'}
            FLAG_COMPLETE = true;

            if nargin > 3
                phi = varargin{1, 2};
                phi_c = varargin{1, 3};

                if phi < 1
                    listSpecies = obj.listSpeciesLean;
                elseif phi >= 1 && phi < phi_c
                    listSpecies = obj.listSpeciesRich;
                else
                    listSpecies = obj.listSpeciesSoot;
                end

            else
                listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb'};
            end

        case 'HC/O2/N2 EXTENDED'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'C2', ...
                         'CH', 'CH3', 'CH4', 'CN', 'H', 'HCN', 'HCO', 'HO2', 'N', 'N2O', ...
                         'NH2', 'NH3', 'NO', 'NO2', 'O', 'OH'};

        case 'HC/O2/N2'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar'};

        case 'HC/O2/N2 RICH'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', ...
                         'C2H4', 'CH', 'CH3', 'CH4', 'CN', 'H', 'HCN', 'HCO', ...
                         'N', 'NH', 'NH2', 'NH3', 'NO', 'O', 'OH'};

        case 'SOOT FORMATION'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
                         'C2', 'C2H4', 'CH', 'CH3', 'CH4', 'CN', 'H', ...
                         'HCN', 'HCO', 'N', 'NH', 'NH2', 'NH3', 'NO', 'O', 'OH'};

        case 'SOOT FORMATION EXTENDED'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
                         'C2', 'C2H', 'C2H2_acetylene', 'C2H2_vinylidene', ...
                         'C2H3_vinyl', 'C2H4', 'C2H5', 'C2H5OH', 'C2H6', ...
                         'C2N2', 'C2O', 'C3', 'C3H3_1_propynl', ...
                         'C3H3_2_propynl', 'C3H4_allene', 'C3H4_propyne', ...
                         'C3H5_allyl', 'C3H6O_acetone', 'C3H6_propylene', ...
                         'C3H8', 'C4', 'C4H2_butadiyne', 'C5', 'C6H2', 'C6H6', ...
                         'C8H18_isooctane', 'CH', 'CH2', 'CH2CO_ketene', ...
                         'CH2OH', 'CH3', 'CH3CHO_ethanal', 'CH3CN', ...
                         'CH3COOH', 'CH3O', 'CH3OH', 'CH4', 'CN', 'COOH', 'H', ...
                         'H2O2', 'HCCO', 'HCHO_formaldehy', 'HCN', 'HCO', ...
                         'HCOOH', 'HNC', 'HNCO', 'HNO', 'HO2', 'N', 'N2O', ...
                         'NCO', 'NH', 'NH2', 'NH2OH', 'NH3', 'NO', 'NO2', ...
                         'O', 'OCCN', 'OH', 'C3O2', 'C4N2', 'CH3CO_acetyl', ...
                         'C4H6_butadiene', 'C4H6_1butyne', 'C4H6_2butyne', ...
                         'C2H4O_ethylen_o', 'CH3OCH3', 'C4H8_1_butene', ...
                         'C4H8_cis2_buten', 'C4H8_isobutene', ...
                         'C4H8_tr2_butene', 'C4H9_i_butyl', 'C4H9_n_butyl', ...
                         'C4H9_s_butyl', 'C4H9_t_butyl', 'C6H5OH_phenol', ...
                         'C6H5O_phenoxy', 'C6H5_phenyl', 'C7H7_benzyl', ...
                         'C7H8', 'C8H8_styrene', 'C10H8_naphthale'};

        case {'AIR', 'DISSOCIATED AIR'}
            listSpecies = {'CO2', 'CO', 'O2', 'N2', 'Ar', 'O', 'O3', ...
                         'N', 'NO', 'NO2', 'NO3', 'N2O', 'N2O3', ...
                         'N2O4', 'N3', 'C'};

        case {'AIR_IONS', 'AIR IONS'}
            listSpecies = {'eminus', 'Ar', 'Arplus', 'C', 'Cplus', 'Cminus', ...
                         'CN', 'CNplus', 'CNminus', 'CNN', 'CO', 'COplus', ...
                         'CO2', 'CO2plus', 'C2', 'C2plus', 'C2minus', 'CCN', ...
                         'CNC', 'OCCN', 'C2N2', 'C2O', 'C3', 'C3O2', 'N', ...
                         'Nplus', 'Nminus', 'NCO', 'NO', 'NOplus', 'NO2', ...
                         'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
                         'N2minus', 'NCN', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
                         'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
                         'O2minus', 'O3'};

        case {'IDEAL_AIR', 'AIR_IDEAL'}
            listSpecies = {'O2', 'N2', 'O', 'O3', 'N', 'NO', 'NO2', 'NO3', 'N2O', ...
                         'N2O3', 'N2O4', 'N3'};

        case 'HYDROGEN'
            listSpecies = {'H2O', 'H2', 'O2', 'N2', 'Ar', 'H', 'HNO', ...
                         'HNO3', 'NH', 'NH2OH', 'NO3', 'N2H2', 'N2O3', 'N3', 'OH', ...
                         'HNO2', 'N', 'NH3', 'NO2', 'N2O', 'N2H4', 'N2O5', 'O', 'O3', ...
                         'HO2', 'NH2', 'H2O2', 'N3H', 'NH2NO2'};

        case {'HYDROGEN_IONS', 'HYDROGEN IONS'}
            listSpecies = {'H2O', 'H2', 'O2', 'N2', 'H', 'OH', 'H2O2', 'H2Oplus', ...
                         'H2minus', 'H2plus', 'H3Oplus', 'HNO', 'HNO2', 'HNO3', 'HO2', ...
                         'HO2minus', 'Hminus', 'Hplus', 'N', 'N2H2', 'N2H4', 'N2O', 'N2O3', ...
                         'N2O5', 'N2Oplus', 'N2minus', 'N2plus', 'N3', 'N3H', 'NH', 'NH2', ...
                         'NH2NO2', 'NH2OH', 'NH3', 'NO2', 'NO2minus', 'NO3', 'NO3minus', ...
                         'NOplus', 'Nminus', 'Nplus', 'O', 'O2minus', 'O2plus', 'O3', ...
                         'Ominus', 'Oplus', 'eminus'};

        case {'HYDROGEN_L', 'HYDROGEN (L)'}
            listSpecies = {'H2O', 'H2', 'O2', 'H', 'OH', 'O', 'O3', 'HO2', ...
                         'H2O2', 'H2bLb', 'O2bLb'};

        case 'HC/O2/N2 PROPELLANTS'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
                         'C2', 'C2H', 'C2H2_acetylene', 'C2H2_vinylidene', ...
                         'C2H3_vinyl', 'C2H4', 'C2H5', 'C2H5OH', 'C2H6', ...
                         'C2N2', 'C2O', 'C3', 'C3H3_1_propynl', ...
                         'C3H3_2_propynl', 'C3H4_allene', 'C3H4_propyne', ...
                         'C3H5_allyl', 'C3H6O_acetone', 'C3H6_propylene', ...
                         'C3H8', 'C4', 'C4H2_butadiyne', 'C5', 'C6H2', 'C6H6', ...
                         'C8H18_isooctane', 'CH', 'CH2', 'CH2CO_ketene', ...
                         'CH2OH', 'CH3', 'CH3CHO_ethanal', 'CH3CN', ...
                         'CH3COOH', 'CH3O', 'CH3OH', 'CH4', 'CN', 'COOH', 'H', ...
                         'H2O2', 'HCCO', 'HCHO_formaldehy', 'HCN', 'HCO', ...
                         'HCOOH', 'HNC', 'HNCO', 'HNO', 'HO2', 'N', 'N2O', ...
                         'NCO', 'NH', 'NH2', 'NH2OH', 'NH3', 'NO', 'NO2', ...
                         'O', 'OCCN', 'OH', 'C3O2', 'C4N2', 'RP_1', 'H2bLb', ...
                         'O2bLb'};

        case 'SI/HC/O2/N2 PROPELLANTS'
            listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
                         'C2', 'C2H4', 'CH', 'CH3', 'CH4', 'CN', 'H', ...
                         'H2O2', 'HCN', 'HCO', 'N', 'NH', 'NH2', 'NH3', 'NO', 'O', 'OH', ...
                         'O2bLb', 'Si', 'SiH', 'SiH2', 'SiH3', 'SiH4', 'SiO2', 'SiO', ...
                         'SibLb', 'SiO2bLb', 'Si2'};
    end

end

function listSpeciesFormula = getFormula(listSpecies, database)
    % Get chemical formula from the database (DB)
    for i = length(listSpecies):-1:1
        listSpeciesFormula{i} = database.species.(listSpecies{i}).formula;
    end

end
