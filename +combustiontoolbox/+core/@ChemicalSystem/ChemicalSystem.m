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
        FLAG_INITIALIZE = true
    end

    methods
        
        [LS, ind_elements_DB] = findProducts(obj, listSpecies, varargin)
        [obj, LS, listSpeciesFormula] = setListSpecies(obj, database, varargin)
        
        function obj = ChemicalSystem(database, varargin)
            % Constructor

            % Definitions
            defaultListSpecies = [];
            
            % Check additional inputs
            if nargin > 1
                varargin = [{'listSpecies'}, varargin(:)'];
            end

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'database'); % , @(x) isa(x, 'combustiontoolbox.databases.NasaDatabase') || isa(x, 'combustiontoolbox.databases.BurcatDatabase')
            addParameter(ip, 'listSpecies', defaultListSpecies); % , @(x) ischar(x) || iscell(x)
            addParameter(ip, 'FLAG_BURCAT', obj.FLAG_BURCAT)
            addParameter(ip, 'FLAG_ION', obj.FLAG_ION)
            parse(ip, database, varargin{:});
            
            % Assign properties
            obj.database = database;
            obj.listSpecies = ip.Results.listSpecies;
            obj.FLAG_BURCAT = ip.Results.FLAG_BURCAT;
            obj.FLAG_ION = ip.Results.FLAG_ION;
            
            % Check if listSpecies is defined
            if isempty(obj.listSpecies)
                obj.FLAG_INITIALIZE = false;
                return
            end
            
            % Initialize chemical system
            obj = initialization(obj);
        end

        function obj = plus(obj, obj2)
            % Overload the plus operator to add propertiesMatrix of two ChemicalSystem objects
            obj.propertiesMatrix(:, obj.ind_ni:obj.numProperties) = obj.propertiesMatrix(:, obj.ind_ni:obj.numProperties) + obj2.propertiesMatrix(:, obj.ind_ni:obj.numProperties);
        end

        function obj = initialization(obj)
            % Initialize chemical system

            % Set list species
            [obj, ~, obj.listSpeciesFormula] = obj.setListSpecies(obj.database, obj.listSpecies);
            
            % Set species
            obj = obj.get_species(obj.database);

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

            % Check if chemical system is initialized
            if ~obj.FLAG_INITIALIZE

                for i = 1:length(species)

                    if ~strcmp(obj.listSpecies, species(i))
                        obj.species.(species{i}) = obj.database.species.(species{i});
                        obj.listSpecies = [obj.listSpecies, species(i)];
                    end
    
                end

                % Include all possible combinations in the database as
                % products
                obj.listSpecies = unique([obj.listSpecies, findProducts(obj, obj.listSpecies)], 'stable');
                
                % Initialize chemical system
                initialization(obj);
                return
            end

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
            %     obj (ChemicalSystem): 
            %
            % Examples:
            %     obj = obj.set_propertiesMatrix({'N2', 'O2'}, [3.76, 1], 300)
            %     obj = obj.set_propertiesMatrix({'N2', 'O2'}, [3.76, 1], 300, [1, 2])
            
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
            
            % Definitions
            if isempty(obj.propertiesMatrix)
                M0 = zeros(obj.numSpecies, obj.numProperties);
            else
                M0 = obj.propertiesMatrix;
            end

            % Get index species
            index = obj.indexSpecies;

            % Fill properties matrix
            for i = obj.numSpecies:-1:1
                M0(index(i), obj.ind_W) = obj.species.(obj.listSpecies{i}).W; % [g/mol]
                M0(index(i), obj.ind_hfi) = obj.species.(obj.listSpecies{i}).hf; % [J/mol]
                M0(index(i), obj.ind_efi) = obj.species.(obj.listSpecies{i}).ef; % [J/mol]
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
                propertiesMatrix(index(i), obj.ind_hi) = species_h0(species{i}, T, obj.species); % [J/mol]
                propertiesMatrix(index(i), obj.ind_cpi) = species_cP(species{i}, T, obj.species); % [J/mol-K]
                propertiesMatrix(index(i), obj.ind_si) = species_s0(species{i}, T, obj.species); % [J/mol-K]
            end

        end
            
        function propertiesMatrix = fillPropertiesMatrixFast(obj, propertiesMatrix, species, moles, T, index)
            % Fill properties matrix
            propertiesMatrix(index, obj.ind_ni) = moles(index); % [mol]
            
            % for i = length(index):-1:1
            %     propertiesMatrix(index(i), obj.ind_hi) = obj.species.(species{i}).get_h0(T); % [J/mol]
            %     propertiesMatrix(index(i), obj.ind_cpi) = obj.species.(species{i}).get_cp(T); % [J/mol-K]
            %     propertiesMatrix(index(i), obj.ind_si) = obj.species.(species{i}).get_s0(T); % [J/mol-K]
            % end

            for i = length(index):-1:1
                propertiesMatrix(index(i), obj.ind_hi) = species_h0(species{i}, T, obj.species); % [J/mol]
                propertiesMatrix(index(i), obj.ind_cpi) = species_cP(species{i}, T, obj.species); % [J/mol-K]
                propertiesMatrix(index(i), obj.ind_si) = species_s0(species{i}, T, obj.species); % [J/mol-K]
            end

        end

        function propertiesMatrix = fillPropertiesMatrixFastH0(obj, propertiesMatrix, species, moles, T, index, h0)
            % Fill properties matrix
            propertiesMatrix(index, obj.ind_ni) = moles(index); % [mol]
            propertiesMatrix(index, obj.ind_hi) = h0(index); % [J/mol]

            for i = length(index):-1:1
                propertiesMatrix(index(i), obj.ind_cpi) = species_cP(species{i}, T, obj.species); % [J/mol-K]
                propertiesMatrix(index(i), obj.ind_si) = species_s0(species{i}, T, obj.species); % [J/mol-K]
            end

        end

    end

end