classdef ChemicalSystem < handle & matlab.mixin.Copyable

    properties
        species                % Struct with class Species
        listSpecies            % List of species
        listElements           % List of elements
        stoichiometricMatrix   % Stoichiometric matrix
        propertiesMatrix       % Properties matrix
        molesPhaseMatrix       % Matrix [moles, phase]
        % * Index values
        indexGas               % Indeces gaseous species
        indexCondensed         % Indeces condensed species
        indexCryogenic         % Indeces cryogenic liquified species
        indexIons              % Indeces ionized species in species
        indexReact             % Indeces react species
        indexFrozen            % Indeces inert/frozen species
        % * List of species for a complete combustion
        listSpeciesLean = {'CO2', 'H2O', 'N2', 'Ar', 'O2'}                % List of species for a lean complete combustion (equivalence ratio < 1)
        listSpeciesRich = {'CO2', 'H2O', 'N2', 'Ar', 'CO', 'H2'}          % List of species for a rich complete combustion (equivalence ratio > 1)
        listSpeciesSoot = {'N2', 'Ar', 'CO', 'H2', 'Cbgrb', 'CO2', 'H2O'} % List of species for a roch complete combustion with soot formation  (equivalence ratio > equivalence ratio soot)
        % * Flags
        FLAG_COMPLETE = false % Flag indicating to compute chemical equilibrium considering a complete combustion
        FLAG_BURCAT = false   % Find all the combinations of species from the database (without BURCAT's DB) that can appear as products for the given list of reactants
        FLAG_ION = false      % Flag indicating to include ionized species in the automatic finder of species
        FLAG_CONDENSED = true % Flag indicating to include condensed species  
    end

    properties (Access = public, Constant = true)
        defaultListSpecies = 'Soot formation'
        ind_hfi   = 1 % Index enthalpy of formation
        ind_efi   = 2 % Index internal energy of formation
        ind_W     = 3 % Index molecular weight
        ind_phase = 4 % Index phase
        ind_ni    = 5 % Index number of moles
        ind_hi    = 6 % Index enthalpy
        ind_cpi   = 7 % Index specific heat at constant pressure
        ind_si    = 8 % Index entropy
        numProperties = 8 % Number of properties in propertiesMatrix
    end

    properties (Dependent)
        numSpecies
        numSpeciesGas
        numElements
        indexElements
        indexSpecies
        ind_C
        ind_H
        ind_O
        ind_N
        ind_E % Index electron
        ind_S % Index 
        ind_Si % Index
    end

    properties (Hidden)
        propertiesMatrixFuel
        propertiesMatrixOxidizer
        listProducts
        indexProducts
        oxidizerReferenceIndex
        oxidizerReferenceAtomsO
    end
    
    properties (Access = private, Hidden)
        listSpeciesFormula
        database
    end

    methods
        
        [LS, ind_elements_DB] = findProducts(obj, listSpecies, varargin)

        function obj = ChemicalSystem(database, listSpecies, varargin)

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'database'); % , @(x) isa(x, 'combustiontoolbox.databases.NasaDatabase') || isa(x, 'combustiontoolbox.databases.BurcatDatabase')
            addRequired(ip, 'listSpecies'); % , @(x) ischar(x) || iscell(x)
            parse(ip, database, listSpecies, varargin{:});
            
            % Assign memory reference of the database
            obj.database = database;

            % Set list species
            [obj, ~, obj.listSpeciesFormula] = obj.list_species(database, ip.Results.listSpecies);
            
            % Set species
            obj = obj.get_species(database);

            % Set contained elements
            obj = obj.containedElements();
            
            % Sort species: first gaseous species, secondly condensed species
            obj = obj.list_phase_species();

            % Set stoichiometric matrix
            obj = obj.setStoichiometricMatrix();

            % Initialize the properties matrix
            obj = obj.setPropertiesMatrixInitialize();

            % Assign listSpecies with the species to be considered in the
            % chemical transformation
            obj.listProducts = obj.listSpecies;
        end
        
        function system = getSystemProducts(obj)
            % Set List of Species to List of Products
            
            % Check if all species are considered as possible products
            FLAG_SAME = all(ismember(obj.indexSpecies, obj.indexProducts));
            
            % Copy memory reference to the chemical system
            if FLAG_SAME
                system = obj;
                return
            end

            % Copy chemical system
            system = obj.copy();

            % Remove ionized species if TP is below T_ions
            % if any(system.indexIons) && temperature < temperatureIons
            %     system.indexListProducts(system.indexIons) = [];
            %     system.indexElements = [];
            % end
            
            % Initialization
            system.indexGas = []; system.indexCondensed = [];
            system.indexCryogenic = []; system.indexIons = [];
            % Set list of species for calculations
            system.listSpecies = system.listSpecies(system.indexProducts);
            % Establish cataloged list of species according to the state of the phase
            system = system.list_phase_species();
            % Update stoichiometric matrix
            system.stoichiometricMatrix = system.stoichiometricMatrix(system.indexProducts, :);
            % Update property matrix
            system.propertiesMatrix = system.propertiesMatrix(system.indexProducts, :);
            % Update compostion matrix
            system.molesPhaseMatrix = system.molesPhaseMatrix(system.indexProducts, :);
        end

        function obj = checkSpecies(obj, species)
            %
            
            % Initialization
            FLAG_ADDED_SPECIES = false;

            % Check if species is defined in the chemical system
            for i = 1:length(species)

                if ~strcmp(obj.listSpecies, species(i))
                    obj.species.(species{i}) = obj.database.species.(species{i});
                    obj.listSpecies = [obj.listSpecies, species(i)];
                    obj.listSpeciesFormula = [obj.listSpeciesFormula, obj.species.(species{i}).formula];
                    FLAG_ADDED_SPECIES = true;
                end

            end

            if ~FLAG_ADDED_SPECIES
                return
            end

            % Set contained elements
            obj = obj.containedElements();
            
            % Sort species: first gaseous species, secondly condensed species
            obj = obj.list_phase_species();

            % Set stoichiometric matrix
            obj = obj.setStoichiometricMatrix();

            % Initialize the properties matrix
            obj = obj.setPropertiesMatrixInitialize();
        end

        function value = get.numSpecies(obj)
            value = length(obj.listSpecies);
        end

        function value = get.numSpeciesGas(obj)
            value = length(obj.indexGas);
        end

        function value = get.numElements(obj)
            value = length(obj.listElements);
        end

        function value = get.ind_C(obj)
            value = find(ismember(obj.listElements, 'C'));
        end

        function value = get.ind_H(obj)
            value = find(ismember(obj.listElements, 'H'));
        end

        function value = get.ind_O(obj)
            value = find(ismember(obj.listElements, 'O'));
        end

        function value = get.ind_N(obj)
            value = find(ismember(obj.listElements, 'N'));
        end

        function value = get.ind_E(obj)
            value = find(ismember(obj.listElements, 'E'));
        end

        function value = get.ind_S(obj)
            value = find(ismember(obj.listElements, 'S'));
        end

        function value = get.ind_Si(obj)
            value = find(ismember(obj.listElements, 'SI'));
        end
    
        function value = get.indexSpecies(obj)
            value = [obj.indexGas, obj.indexCondensed];
        end

        function value = get.indexElements(obj)
            value = 1:1:obj.numElements;
        end

        function obj = get_species(obj, database)
            
            for i = obj.numSpecies:-1:1
                obj.species.(obj.listSpecies{i}) = database.species.(obj.listSpecies{i});
            end

        end

        function obj = set_propertiesMatrix(obj, species, moles, T, varargin)
            % Fill the properties matrix with the data of the mixture
            %
            % Args:
            %     obj (ChemicalSystem): 
            %     species (cell): Species contained in the system
            %     moles (float): Moles of the species in the mixture [mol]
            %     T (float): Temperature [K]
            %
            % Optional Args:
            %     ind (float): Vector with the indexes of the species to fill the properties matrix   
            %
            % Returns:
            %     properties_matrix (float): Properties matrix
            %
            % Examples:
            %     obj = obj.set_species({'N2', 'O2'}, [3.76, 1], 300)
            %     obj = obj.set_species({'N2', 'O2'}, [3.76, 1], 300, [1, 2])
            
            % Fill properties matrix
            if nargin < 5
                obj.propertiesMatrix = obj.fillPropertiesMatrix(obj, obj.propertiesMatrix, species, moles, T);
                return
            elseif nargin < 6
                index = varargin{1};
                obj.propertiesMatrix = obj.fillPropertiesMatrixFast(obj, obj.propertiesMatrix, species(index), moles, T, index);
                return
            end
    
            index = varargin{1};
            h0 = varargin{2};
            obj.propertiesMatrix = obj.fillPropertiesMatrixFastH0(obj, obj.propertiesMatrix, species(index), moles, T, index, h0);
        end

        function obj = clean(obj)
            % Set temperature-dependent matrix properties to zero
            obj.propertiesMatrix(:, 5:end) = 0;
        end

        function obj = cleanMoles(obj)
            % Set temperature-dependent matrix properties to zero
            obj.propertiesMatrix(:, obj.ind_ni) = 0;
        end

        function [obj, LS, listSpeciesFormula] = list_species(obj, database, varargin)
            % Set list of species in the mixture (products)
            %
            % Predefined list of species:
            %     * SOOT FORMATION (default)
            %     * COMPLETE
            %     * HC/O2/N2 EXTENDED
            %     * SOOT FORMATION EXTENDED
            %     * NASA ALL
            %     * NASA ALL CONDENSED
            %     * NASA ALL IONS
            %     * AIR, DISSOCIATED AIR
            %     * AIR IONS, AIR_IONS
            %     * IDEAL_AIR, AIR_IDEAL
            %     * HYDROGEN
            %     * HYDROGEN_L, HYDROGEN (L)
            %     * HC/O2/N2 PROPELLANTS
            %     * SI/HC/O2/N2 PROPELLANTS
            %
            % Optional Args:
            %     * self (struct): Data of the mixture, conditions, and databases
            %     * LS (cell): Name list species / list of species
            %     * phi (float): Equivalence ratio
            %     * phi_c (float): Equivalence ratio in which theoretically appears soot
            %
            % Returns:
            %     Tuple containing
            %
            %     * self (struct): Data of the mixture, conditions, and databases
            %     * LS (cell): List of species
            %
            % Examples:
            %     * LS = list_species('soot formation');
            %     * [self, LS] = list_species(self, 'soot formation');
            %     * [self, LS] = list_species(self, 'complete', 1.5, 2.5);
        
            % Unpack inputs
            [obj, LS, FLAG] = unpack(obj, varargin{:});

            % Set ListSpecies (LS)
            if ~isempty(LS)
        
                switch upper(LS)
                    case {'COMPLETE', 'COMPLETE REACTION'}
                        obj.FLAG_COMPLETE = true;
        
                        if nargin > 2
                            phi = varargin{1, 3};
                            phi_c = varargin{1, 4};
        
                            if phi < 1
                                obj.listSpecies = obj.listSpeciesLean;
                            elseif phi >= 1 && phi < phi_c
                                obj.listSpecies = obj.listSpeciesRich;
                            else
                                obj.listSpecies = obj.listSpeciesSoot;
                            end
        
                        else
                            obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb'};
                        end
        
                    case 'HC/O2/N2 EXTENDED'
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'C2', ...
                                     'CH', 'CH3', 'CH4', 'CN', 'H', 'HCN', 'HCO', 'HO2', 'N', 'N2O', ...
                                     'NH2', 'NH3', 'NO', 'NO2', 'O', 'OH'};
        
                    case 'HC/O2/N2'
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar'};
        
                    case 'HC/O2/N2 RICH'
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', ...
                                     'C2H4', 'CH', 'CH3', 'CH4', 'CN', 'H', 'HCN', 'HCO', ...
                                     'N', 'NH', 'NH2', 'NH3', 'NO', 'O', 'OH'};
        
                    case 'SOOT FORMATION'
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
                                     'C2', 'C2H4', 'CH', 'CH3', 'CH4', 'CN', 'H', ...
                                     'HCN', 'HCO', 'N', 'NH', 'NH2', 'NH3', 'NO', 'O', 'OH'};
        
                    case 'SOOT FORMATION EXTENDED'
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
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
        
                    case 'NASA ALL'
                        obj.listSpecies = find_species_LS(fieldnames(database.species), ...
                            {'C', 'N', 'O', 'minus', 'plus', 'Ar', 'H'}, 'any', ...
                            {'I', 'S', 'L', 'T', 'P', 'F', 'ab', 'W', 'AL','He' ...
                                'Z', 'X', 'R', 'Os', 'Cr', 'Br', 'G', 'K','Li' ...
                                'U', 'Co', 'Cu', 'B', 'V', 'Ni', 'Na', 'Mg','Hg' ...
                                'Mo', 'Ag', 'Nb', 'Cb', 'CL', 'D', 'T', 'M','minus' ...
                                'Ca', 'Cs', 'Ne', 'Cd', 'Mn', 'cr', 'plus'}, 'all');
        
                    case {'NASA ALL CONDENSED', 'NASA ALL_CONDENSED', 'NASA_ALL_CONDENSED'}
                        obj.listSpecies = find_species_LS(fieldnames(database.species), ...
                            {'C', 'N', 'O', 'Ar', 'H'}, 'any', ...
                            {'I', 'S', 'T', 'P', 'F', 'ab', 'W', 'AL','He' ...
                                'Z', 'X', 'R', 'Os', 'Cr', 'Br', 'G', 'K','Li' ...
                                'U', 'Co', 'Cu', 'B', 'V', 'Ni', 'Na', 'Mg','Hg' ...
                                'Mo', 'Ag', 'Nb', 'Cb', 'CL', 'D', 'T', 'M', 'minus', ...
                                'Ca', 'Cs', 'Ne', 'Cd', 'Mn', 'plus'}, 'all');
        
                    case {'NASA ALL IONS', 'NASA ALL_IONS', 'NASA_ALL_IONS'}
                        obj.listSpecies = find_species_LS(fieldnames(database.species), ...
                            {'C', 'N', 'O', 'minus', 'plus', 'Ar', 'H'}, 'any', ...
                            {'I', 'S', 'L', 'T', 'P', 'F', 'ab', 'W', 'AL','He' ...
                                'Z', 'X', 'R', 'Os', 'Cr', 'Br', 'G', 'K','Li' ...
                                'U', 'Co', 'Cu', 'B', 'V', 'Ni', 'Na', 'Mg','Hg' ...
                                'Mo', 'Ag', 'Nb', 'Cb', 'CL', 'D', 'T', 'M', ...
                                'Ca', 'Cs', 'Ne', 'Cd', 'Mn', 'cr', 'cyclo'}, 'all');
        
                    case {'AIR', 'DISSOCIATED AIR'}
                        obj.listSpecies = {'CO2', 'CO', 'O2', 'N2', 'Ar', 'O', 'O3', ...
                                     'N', 'NO', 'NO2', 'NO3', 'N2O', 'N2O3', ...
                                     'N2O4', 'N3', 'C'};
        
                    case {'AIR_IONS', 'AIR IONS'}
                        obj.listSpecies = {'eminus', 'Ar', 'Arplus', 'C', 'Cplus', 'Cminus', ...
                                     'CN', 'CNplus', 'CNminus', 'CNN', 'CO', 'COplus', ...
                                     'CO2', 'CO2plus', 'C2', 'C2plus', 'C2minus', 'CCN', ...
                                     'CNC', 'OCCN', 'C2N2', 'C2O', 'C3', 'C3O2', 'N', ...
                                     'Nplus', 'Nminus', 'NCO', 'NO', 'NOplus', 'NO2', ...
                                     'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
                                     'N2minus', 'NCN', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
                                     'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
                                     'O2minus', 'O3'};
        
                    case {'IDEAL_AIR', 'AIR_IDEAL'}
                        obj.listSpecies = {'O2', 'N2', 'O', 'O3', 'N', 'NO', 'NO2', 'NO3', 'N2O', ...
                                     'N2O3', 'N2O4', 'N3'};
        
                    case 'HYDROGEN'
                        obj.listSpecies = {'H2O', 'H2', 'O2', 'N2', 'Ar', 'H', 'HNO', ...
                                     'HNO3', 'NH', 'NH2OH', 'NO3', 'N2H2', 'N2O3', 'N3', 'OH', ...
                                     'HNO2', 'N', 'NH3', 'NO2', 'N2O', 'N2H4', 'N2O5', 'O', 'O3', ...
                                     'HO2', 'NH2', 'H2O2', 'N3H', 'NH2NO2'};
        
                    case {'HYDROGEN_IONS', 'HYDROGEN IONS'}
                        obj.listSpecies = {'H2O', 'H2', 'O2', 'N2', 'H', 'OH', 'H2O2', 'H2Oplus', ...
                                     'H2minus', 'H2plus', 'H3Oplus', 'HNO', 'HNO2', 'HNO3', 'HO2', ...
                                     'HO2minus', 'Hminus', 'Hplus', 'N', 'N2H2', 'N2H4', 'N2O', 'N2O3', ...
                                     'N2O5', 'N2Oplus', 'N2minus', 'N2plus', 'N3', 'N3H', 'NH', 'NH2', ...
                                     'NH2NO2', 'NH2OH', 'NH3', 'NO2', 'NO2minus', 'NO3', 'NO3minus', ...
                                     'NOplus', 'Nminus', 'Nplus', 'O', 'O2minus', 'O2plus', 'O3', ...
                                     'Ominus', 'Oplus', 'eminus'};
        
                    case {'HYDROGEN_L', 'HYDROGEN (L)'}
                        obj.listSpecies = {'H2O', 'H2', 'O2', 'H', 'OH', 'O', 'O3', 'HO2', ...
                                     'H2O2', 'H2bLb', 'O2bLb'};
        
                    case 'HC/O2/N2 PROPELLANTS'
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
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
                        obj.listSpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar', 'Cbgrb', ...
                                     'C2', 'C2H4', 'CH', 'CH3', 'CH4', 'CN', 'H', ...
                                     'H2O2', 'HCN', 'HCO', 'N', 'NH', 'NH2', 'NH3', 'NO', 'O', 'OH', ...
                                     'O2bLb', 'Si', 'SiH', 'SiH2', 'SiH3', 'SiH4', 'SiO2', 'SiO', ...
                                     'SibLb', 'SiO2bLb', 'Si2'};
                    otherwise
        
                        if iscell(LS)
                            obj.listSpecies = LS;
                        else
                            obj.listSpecies = {LS};
                        end
        
                end
        
            end
        
            obj.listSpecies = unique(obj.listSpecies, 'stable');
            LS = obj.listSpecies;
        
            if FLAG
                self = obj.listSpecies;
                return
            end
        
            listSpeciesFormula = get_formula(obj.listSpecies, database);
        
            if any(obj.getIndexIons(obj.listSpecies))
                obj.FLAG_ION = true;
            end

            % SUB-PASS FUNCTIONS
            function listSpeciesFormula = get_formula(listSpecies, database)
                % Get chemical formula from the database (DB)
                for i = length(listSpecies):-1:1
                    listSpeciesFormula{i} = database.species.(listSpecies{i}).formula;
                end
            
            end
            
            function [obj, LS] = unpack_LS(obj, variable)
                % Unpack list of species (LS)
                LS = [];
            
                if iscell(variable)
                    obj.listSpecies = variable;
                else
                    LS = variable;
                end
            
            end
            
            function [obj, LS, FLAG] = unpack(obj, varargin)
                
                % Unpack
                if nargin < 2
                    FLAG = true; % Return variable "LS"
            
                    if ~iscell(varargin{1}) && ~ischar(varargin{1})
                        FLAG = false; % Return variable "self"
                        obj = varargin{1};
                        return
                    end
            
                    [obj, LS] = unpack_LS(obj, varargin{1});
            
                else
                    FLAG = false; % Return variable "self"
                    [obj, LS] = unpack_LS(obj, varargin{1});
                end
            
            end
            
        end

        function obj = check_complete_reaction(obj, equivalenceRatio, equivalenceRatioSoot)
            % Check if the list of species corresponds to "complete_reaction"
            % If FLAG_COMPLETE is true, establish the list of species based on the
            % given equivalence ratio (phi)
            if ~obj.FLAG_COMPLETE
                return
            end
        
            if equivalenceRatio < 1
                LS = obj.LS_lean;
            elseif equivalenceRatio >= 1 && equivalenceRatio < equivalenceRatioSoot
                LS = obj.LS_rich;
            else
                LS = obj.LS_soot;
            end
        
            obj.indexProducts = findIndex(obj.listSpecies, LS);
            obj = obj.reorganizeIndexPhaseSpecies();
        end

        function obj = list_phase_species(obj)
            % Establish cataloged list of species according to the state of the
            % phase (gaseous or condensed). It also obtains the indices of
            % cryogenic liquid species, i.e., liquified gases.
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %     LS (cell):     List of species
            %
            % Returns:
            %     self (struct): Data of the mixture, conditions, and databases
        
            obj = obj.getIndexPhaseSpecies();
            obj.listSpecies = obj.listSpecies([obj.indexGas, obj.indexCondensed]);
            % Reorginize index of gaseous, condensed and cryogenic species
            obj = obj.reorganizeIndexPhaseSpecies();
        end

        function obj = set_react_index(obj, speciesFrozen)
            % Set index of react (non-frozen) and frozen species
            %
            % Args:
            %     obj (ChemicalSystem): 
            %     species (char): Frozen species
            %
            % Returns:
            %     self (struct): Data of the mixture, conditions, and databases
            
            % Import packages
            import combustiontoolbox.utils.findIndex

            % Initialization
            obj.indexReact = 1:obj.numSpecies;
        
            % All species react
            if isempty(speciesFrozen)
                obj.indexFrozen = [];
                return
            end
            
            % Get index frozen species
            index = findIndex(obj.listSpecies, speciesFrozen);
            % Set index frozen species
            obj.indexFrozen = index;

            % Get length initial species
            % N_reactants = length([self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert]);
            % Get length frozen species
            N_frozen = length(index);
            % Check if all the species of the reactants are frozen
            % if N_frozen == N_reactants
            %     self.PD.FLAG_FROZEN = true;
            %     return
            % end
            
            for i = 1:N_frozen
                obj.indexReact(obj.indexReact == index(i)) = [];
            end
        
        end

        function obj = getOxidizerReference(obj, listOxidizer)
            % Get oxidizer of reference for computations with the equivalence ratio
            %
            % Args:
            %     obj (ChemicalSystem): 
            %     listOxidizer (cell):
            %
            % Returns:
            %     obj (ChemicalSystem):     
            
            % Import packages
            import combustiontoolbox.utils.findIndex

            % Check if there are oxidizers in the mixtures
            if isempty(listOxidizer)
                obj.oxidizerReferenceIndex = [];
                obj.oxidizerReferenceAtomsO = NaN;
                return
            end
        
            % If O2 or O2(L) are included as oxidizers these species will be
            % selected as reference oxidizers in this order. Otherwise, the first
            % oxidizer with oxygen as element will be selected.
            if any(ismember(listOxidizer, 'O2'))
                obj.oxidizerReferenceIndex = findIndex(obj.listSpecies, 'O2');
                obj.oxidizerReferenceAtomsO = 2;
            elseif any(ismember(listOxidizer, 'O2bLb'))
                obj.oxidizerReferenceIndex = findIndex(obj.listSpecies, 'O2bLb');
                obj.oxidizerReferenceAtomsO = 2;
            else
                % Get first oxidizer with oxygen as element
                temp_ind = find(contains(listOxidizer, 'O'), 1);
                species = listOxidizer{temp_ind};
                % Find index of reference oxidizer
                obj.oxidizerReferenceIndex = findIndex(obj.listSpecies, species);
                % Find position oxygen element
                temp_ind_O = find(listOxidizer{temp_ind} == 'O');
                % Get position numbers and letters
                [temp_ind_1, temp_ind_2] = regexp(species, '\w\d*');
                % Find position oxygen element in the temp variable index
                temp_ind = find(temp_ind_1 == temp_ind_O);
                % Set number of elements of oxygen in the reference oxidizer
                obj.oxidizerReferenceAtomsO = sscanf(species(temp_ind_1(temp_ind) + 1:temp_ind_2(temp_ind)), '%f');
            end

        end

    end

    methods (Access = private)
        
        function obj = containedElements(obj)
            % Obtain containted elements from the given set of species (reactants and products)
            %
            % Args:
            %     obj (ChemicalSystem): 
            %
            % Returns:
            %     obj (ChemicalSystem): 
        
            L_formula = [];
        
            for k = obj.numSpecies:-1:1
                L_E1 = []; L_E2 = [];
                formula = obj.listSpeciesFormula{k};
        
                [idx0, idxf] = regexp(formula, "[A-Z]{2,}");
        
                for j = length(idxf):-1:1
                    L_E2{j} = formula(idx0(j):idxf(j));
                    formula(idx0(j):idxf(j)) = ' ';
        
                    if isempty(L_E2)
                        L_E2{j} = [];
                    end
        
                end
        
                [~, idxf] = regexp(formula, "[A-Z]{1}");
        
                for j = length(idxf):-1:1
                    L_E1{j} = formula(idxf(j));
                end
        
                L_formula = [L_formula, L_E1, L_E2];
            end
        
            obj.listElements = unique(L_formula);
        end

        function index = getIndexIons(obj, species)
            % Get index of ions for the given list of species
            %
            % Args:
            %     species (cell): List of species
            %
            % Returns:
            %     index (float): Index of ions
        
            index = (contains(species, 'minus') | contains(species, 'plus')) & ~contains(species, 'cyclominus');
        end
        
        function obj = getIndexPhaseSpecies(obj)
            % Get index of gaseous, condensed and cryogenic species
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %
            % Returns:
            %     self (struct): Data of the mixture, conditions, and databases
        
            % Preallocate arrays
            obj.indexGas = zeros(1, obj.numSpecies);
            obj.indexCondensed = zeros(1, obj.numSpecies);
            obj.indexCryogenic = zeros(1, obj.numSpecies);
            
            % Initialization
            gasCount = 0;
            condensedCount = 0;
            cryogenicCount = 0;
        
            % Get indices
            for index = 1:obj.numSpecies
                species = obj.listSpecies{index};
        
                if ~obj.species.(species).phase
                    gasCount = gasCount + 1;
                    obj.indexGas(gasCount) = index;
                else
                    condensedCount = condensedCount + 1;
                    obj.indexCondensed(condensedCount) = index;
        
                    if ~obj.species.(species).Tintervals
                        cryogenicCount = cryogenicCount + 1;
                        obj.indexCryogenic(cryogenicCount) = index;
                    end
                end
            end
        
            % Trim excess zeros from preallocated arrays
            obj.indexGas = obj.indexGas(1:gasCount);
            obj.indexCondensed = obj.indexCondensed(1:condensedCount);
            obj.indexCryogenic = obj.indexCryogenic(1:cryogenicCount);
        
            % Get index of ions
            obj.indexIons = obj.getIndexIons(obj.listSpecies);
        end        

        function obj = reorganizeIndexPhaseSpecies(obj)
            % Reorginize index of gaseous, condensed and cryogenic species
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %     LS (cell): Name list species / list of species
            %
            % Returns:
            %     self (struct): Data of the mixture, conditions, and databases
        
            obj.indexGas = []; obj.indexCondensed = []; obj.indexCryogenic = [];
            obj = obj.getIndexPhaseSpecies();
        end

        function obj = setStoichiometricMatrix(obj)
            % Set stoichiometric matrix

            % Preallocate the stoichiometric matrix
            A0 = zeros(obj.numSpecies, obj.numElements);
            
            % Set stoichiometric matrix
            for i = 1:obj.numSpecies
                obj.species.(obj.listSpecies{i}).elementMatrix = obj.species.(obj.listSpecies{i}).getElementMatrix(obj.listElements);
                A0(i, obj.species.(obj.listSpecies{i}).elementMatrix(1, :)) = obj.species.(obj.listSpecies{i}).elementMatrix(2, :);
            end

            % Set stoichiometric matrix
            obj.stoichiometricMatrix = A0;
        end

        function obj = setPropertiesMatrixInitialize(obj)
            % Initialize properties matrix
            %
            % Args:
            %     self (struct):  Data of the mixture, conditions, and databases
            %
            % Returns:
            %     self (struct):  Data of the mixture, conditions, and databases
            
            % Preallocation
            M0 = zeros(obj.numSpecies, obj.numProperties);

            % Get index species
            index = obj.indexSpecies;

            % Fill properties matrix
            for i = obj.numSpecies:-1:1
                M0(index(i), obj.ind_W) = obj.species.(obj.listSpecies{i}).W; % [g/mol]
                M0(index(i), obj.ind_hfi) = obj.species.(obj.listSpecies{i}).hf / 1000; % [kJ/mol]
                M0(index(i), obj.ind_efi) = obj.species.(obj.listSpecies{i}).ef / 1000; % [kJ/mol]
                M0(index(i), obj.ind_phase) = obj.species.(obj.listSpecies{i}).phase; % [bool]
            end

            % Initialize molesPhaseMatrix [moles, phase]
            obj.molesPhaseMatrix = M0(:, [obj.ind_ni, obj.ind_phase]);
            
            % Set properties matrix
            obj.propertiesMatrix = M0;
        end

    end

    methods (Access = private, Static)
        
        function propertiesMatrix = fillPropertiesMatrix(obj, propertiesMatrix, species, moles, T)
            %
            
            % Import packages
            import combustiontoolbox.utils.findIndex

            % Get index species
            index = findIndex(obj.listSpecies, species);

            % Fill properties matrix
            propertiesMatrix(index, obj.ind_ni) = moles; % [mol]

            for i = length(moles):-1:1
                propertiesMatrix(index(i), obj.ind_hi) = species_h0(species{i}, T, obj.species); % [kJ/mol]
                propertiesMatrix(index(i), obj.ind_cpi) = species_cP(species{i}, T, obj.species); % [J/mol-K]
                propertiesMatrix(index(i), obj.ind_si) = species_s0(species{i}, T, obj.species); % [kJ/mol-K]
            end

        end
            
        function propertiesMatrix = fillPropertiesMatrixFast(obj, propertiesMatrix, species, moles, T, index)
            % Fill properties matrix
            propertiesMatrix(index, obj.ind_ni) = moles(index); % [mol]
            
            % for i = length(index):-1:1
            %     propertiesMatrix(index(i), obj.ind_hi) = obj.species.(species{i}).get_h0(T); % [kJ/mol]
            %     propertiesMatrix(index(i), obj.ind_cpi) = obj.species.(species{i}).get_cp(T); % [J/mol-K]
            %     propertiesMatrix(index(i), obj.ind_si) = obj.species.(species{i}).get_s0(T); % [kJ/mol-K]
            % end

            for i = length(index):-1:1
                propertiesMatrix(index(i), obj.ind_hi) = species_h0(species{i}, T, obj.species); % [kJ/mol]
                propertiesMatrix(index(i), obj.ind_cpi) = species_cP(species{i}, T, obj.species); % [J/mol-K]
                propertiesMatrix(index(i), obj.ind_si) = species_s0(species{i}, T, obj.species); % [kJ/mol-K]
            end

        end

        function propertiesMatrix = fillPropertiesMatrixFastH0(obj, propertiesMatrix, species, moles, T, index, h0)
            % Fill properties matrix
            propertiesMatrix(index, obj.ind_ni) = moles(index); % [mol]
            propertiesMatrix(index, obj.ind_hi) = h0(index); % [kJ/mol]

            for i = length(index):-1:1
                propertiesMatrix(index(i), obj.ind_cpi) = species_cP(species{i}, T, obj.species); % [J/mol-K]
                propertiesMatrix(index(i), obj.ind_si) = species_s0(species{i}, T, obj.species); % [kJ/mol-K]
            end

        end

    end

end