classdef ChemicalSystem < handle & matlab.mixin.Copyable
    % The :mat:func:`ChemicalSystem` class is used to define a chemical system
    % with a list of species and their thermodynamic properties.
    % 
    % The :mat:func:`ChemicalSystem` object can be initialized as follows: ::
    %
    %       system = ChemicalSystem(DB)
    %       system = ChemicalSystem(DB, {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Ar'})
    %       system = ChemicalSystem(DB, 'soot formation extended')
    %
    % Here ``DB`` represents an instance of the :mat:func:`NasaDatabase` or
    % :mat:func:`BurcatDatabase` class. The first case initializes the chemical
    % system with all the species in the database that may appear as products
    % given the initial mixture. The second case initializes the chemical system
    % with the list of species provided. The third case initialize the chemical
    % system with one of the predefined list of species (see :mat:func:`setListSpecies`)
    % that typically appear in hydrocarbon combustion systems.
    %
    % See also: :mat:func:`Database`, :mat:func:`NasaDatabase`, :mat:func:`findProducts`, :mat:func:`setListSpecies`

    properties
        species                 % Struct with Species objects
        listSpecies             % List of species
        listElements            % List of elements
        stoichiometricMatrix    % Stoichiometric matrix
        molecularWeight         % Molecular weights of species [kg/mol]
        phase                   % Phase vector [-]
        formationEnthalpy       % Enthalpy of formation [J/mol]
        formationInternalEnergy % Internal energy of formation [J/mol]
        temperatureMin          % Minimum species temperature [K]
        temperatureMax          % Maximum species temperature [K]
        indexSpecies            % Index of species
        indexGas                % Indices gaseous species
        indexCondensed          % Indices condensed species
        indexCryogenic          % Indices cryogenic liquified species
        indexIons               % Indices ionized species in species
        indexReact              % Indices react species
        indexFrozen             % Indices inert/frozen species
        listSpeciesLean = {'CO2', 'H2O', 'N2', 'Ar', 'O2'}       % List of species for a lean complete combustion (equivalence ratio < 1)
        listSpeciesRich = {'CO2', 'H2O', 'N2', 'Ar', 'CO', 'H2'} % List of species for a rich complete combustion (equivalence ratio > 1)
        listSpeciesSoot = {'N2', 'Ar', 'CO', 'H2', 'Cbgrb'}      % List of species for a roch complete combustion with soot formation  (equivalence ratio > equivalence ratio soot)
        FLAG_COMPLETE = false   % Flag indicating to compute chemical equilibrium considering a complete combustion
        FLAG_BURCAT = false     % Find all the combinations of species from the database (without BURCAT's DB) that can appear as products for the given list of reactants
        FLAG_ION = false        % Flag indicating to include ionized species in the automatic finder of species
        FLAG_CONDENSED = true   % Flag indicating to include condensed species  
    end

    properties (Dependent)
        numSpecies    % Number of species
        numSpeciesGas % Number of gaseous species
        numElements   % Number of elements
        indexElements % Index of elements
    end

    properties (Hidden)
        database
        listProducts
        oxidizerReferenceIndex
        oxidizerReferenceAtomsO
        ind_C         % Index carbon
        ind_H         % Index hydrogen
        ind_O         % Index oxygen
        ind_N         % Index nitrogen
        ind_E         % Index electron
        ind_S         % Index sulfur
        ind_Si        % Index silicon
        ind_B         % Index boron
    end
    
    properties (Access = private, Hidden)
        listSpeciesFormula     % List of species formulas
        FLAG_INITIALIZE = true % Flag indicating if the chemical system is initialized
    end

    methods
        
        [LS, ind_elements_DB] = findProducts(obj, listSpecies, varargin)
        [obj, listSpecies, listSpeciesFormula] = setListSpecies(obj, varargin)

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
            addOptional(ip, 'listSpecies', defaultListSpecies); % , @(x) ischar(x) || iscell(x)
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

        function obj = initialization(obj)
            % Initialize chemical system

            % Set list species
            setListSpecies(obj, obj.listSpecies);
            
            % Set species
            setSpecies(obj, obj.database);

            % Set contained elements
            setContainedElements(obj);
            
            % Sort species: first gaseous species, secondly condensed species
            sortListSpecies(obj);

            % Set stoichiometric matrix
            setStoichiometricMatrix(obj);

            % Set static species properties
            setStaticSpeciesProperties(obj);

            % Assign listSpecies with the species to be considered in the
            % chemical transformation
            obj.listProducts = obj.listSpecies;
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
            setContainedElements(obj);
            
            % Sort species: first gaseous species, secondly condensed species
            sortListSpecies(obj);

            % Set stoichiometric matrix
            setStoichiometricMatrix(obj);

            % Set static species properties
            setStaticSpeciesProperties(obj);
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
    
        function value = get.indexSpecies(obj)
            value = [obj.indexGas, obj.indexCondensed];
        end

        function value = get.indexElements(obj)
            value = 1:1:obj.numElements;
        end

        function obj = setSpecies(obj, database)
            
            for i = obj.numSpecies:-1:1
                obj.species.(obj.listSpecies{i}) = database.species.(obj.listSpecies{i});
            end

        end

        function thermo = evaluateSpeciesThermo(obj, T, index)
            % Evaluate all temperature-dependent species properties without mutating the system
            %
            % Args:
            %     T (float): Temperature [K]
            %     index (float): Species indices to evaluate
            %
            % Returns:
            %     thermo (struct): Temperature-dependent species properties

            if nargin < 3 || isempty(index)
                index = 1:obj.numSpecies;
            end

            index = index(:).';
            T = reshape(T, 1, []);
            [cpCurves, h0Curves, s0Curves, g0Curves] = combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache(obj, index, 1);
            numSpecies = numel(index);
            numTemperatures = numel(T);
            thermo.index = index;
            thermo.enthalpy = zeros(numSpecies, numTemperatures);
            thermo.heatCapacityPressure = zeros(numSpecies, numTemperatures);
            thermo.entropy = zeros(numSpecies, numTemperatures);
            thermo.gibbs = zeros(numSpecies, numTemperatures);

            for i = 1:numSpecies
                thermo.enthalpy(i, :) = h0Curves{i}(T);
                thermo.heatCapacityPressure(i, :) = cpCurves{i}(T);
                thermo.entropy(i, :) = s0Curves{i}(T);
                thermo.gibbs(i, :) = g0Curves{i}(T);
            end
        end

        function [enthalpy, heatCapacityPressure, entropy] = evaluateSpeciesThermoHCPS(obj, T, index)
            % Evaluate h, cp, and s for the requested species

            if nargin < 3 || isempty(index)
                index = 1:obj.numSpecies;
            end

            index = index(:).';
            T = reshape(T, 1, []);
            [cpCurves, h0Curves, s0Curves] = combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache(obj, index, 2);
            numSpecies = numel(index);
            numTemperatures = numel(T);
            enthalpy = zeros(numSpecies, numTemperatures);
            heatCapacityPressure = zeros(numSpecies, numTemperatures);
            entropy = zeros(numSpecies, numTemperatures);

            for i = 1:numSpecies
                enthalpy(i, :) = h0Curves{i}(T);
                heatCapacityPressure(i, :) = cpCurves{i}(T);
                entropy(i, :) = s0Curves{i}(T);
            end
        end

        function [heatCapacityPressure, entropy] = evaluateSpeciesThermoCPS(obj, T, index)
            % Evaluate cp and s for the requested species

            if nargin < 3 || isempty(index)
                index = 1:obj.numSpecies;
            end

            index = index(:).';
            T = reshape(T, 1, []);
            [cpCurves, ~, s0Curves] = combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache(obj, index, 3);
            numSpecies = numel(index);
            numTemperatures = numel(T);
            heatCapacityPressure = zeros(numSpecies, numTemperatures);
            entropy = zeros(numSpecies, numTemperatures);

            for i = 1:numSpecies
                heatCapacityPressure(i, :) = cpCurves{i}(T);
                entropy(i, :) = s0Curves{i}(T);
            end
        end

        function enthalpy = evaluateSpeciesThermoH(obj, T, index)
            % Evaluate h for the requested species

            if nargin < 3 || isempty(index)
                index = 1:obj.numSpecies;
            end

            index = index(:).';
            T = reshape(T, 1, []);
            [~, h0Curves] = combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache(obj, index, 4);
            numSpecies = numel(index);
            numTemperatures = numel(T);
            enthalpy = zeros(numSpecies, numTemperatures);

            for i = 1:numSpecies
                enthalpy(i, :) = h0Curves{i}(T);
            end
        end

        function [enthalpy, gibbs] = evaluateSpeciesThermoHG(obj, T, index)
            % Evaluate h and g for the requested species

            if nargin < 3 || isempty(index)
                index = 1:obj.numSpecies;
            end

            index = index(:).';
            T = reshape(T, 1, []);
            [~, h0Curves, ~, g0Curves] = combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache(obj, index, 5);
            numSpecies = numel(index);
            numTemperatures = numel(T);
            enthalpy = zeros(numSpecies, numTemperatures);
            gibbs = zeros(numSpecies, numTemperatures);

            for i = 1:numSpecies
                enthalpy(i, :) = h0Curves{i}(T);
                gibbs(i, :) = g0Curves{i}(T);
            end
        end

        function gibbs = evaluateSpeciesThermoG(obj, T, index)
            % Evaluate g for the requested species

            if nargin < 3 || isempty(index)
                index = 1:obj.numSpecies;
            end

            index = index(:).';
            T = reshape(T, 1, []);
            [~, ~, ~, g0Curves] = combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache(obj, index, 6);
            numSpecies = numel(index);
            numTemperatures = numel(T);
            gibbs = zeros(numSpecies, numTemperatures);

            for i = 1:numSpecies
                gibbs(i, :) = g0Curves{i}(T);
            end
        end

        function obj = sortListSpecies(obj)
            % Establish cataloged list of species according to the state of the
            % phase (gaseous or condensed). It also obtains the indices of
            % cryogenic liquid species, i.e., liquified gases.
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the list of species sorted
        
            obj = obj.setIndexPhaseSpecies();
            obj.listSpecies = obj.listSpecies([obj.indexGas, obj.indexCondensed]);

            % Reorginize index of gaseous, condensed and cryogenic species
            obj = obj.sortIndexPhaseSpecies();
        end

        function obj = setReactIndex(obj, speciesFrozen)
            % Set index of react (non-frozen) and frozen species
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %     speciesFrozen (char): Frozen chemical species
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the index of react and frozen species

            % Initialization
            obj.indexReact = 1:obj.numSpecies;
        
            % All species react
            if isempty(speciesFrozen)
                obj.indexFrozen = [];
                return
            end
            
            % Get index frozen species
            index = combustiontoolbox.utils.findIndex(obj.listSpecies, speciesFrozen);

            % Set index frozen species
            obj.indexFrozen = index;

            % Get length frozen species
            N_frozen = length(index);
            
            for i = 1:N_frozen
                obj.indexReact(obj.indexReact == index(i)) = [];
            end
        
        end

        function obj = setOxidizerReference(obj, listOxidizer)
            % Set oxidizer of reference for computations with the equivalence ratio
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %     listOxidizer (cell): List of oxidizers
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the reference oxidizer   
            
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

        function value = isIonized(obj, varargin)
            % Get boolean value indicating if the species is ionized
            %
            % Args:
            %     species (char): Chemical species
            %
            % Returns:
            %     value (bool): Boolean value indicating if the species is ionized
        
            % Definitions
            species = obj.listSpecies;

            % Check additional inputs
            if nargin > 1
                species = varargin{1};
            end

            % Check if species are ionized
            value = (contains(species, 'minus') | contains(species, 'plus')) & ~contains(species, 'cyclominus');
        end

        function charges = getCharges(obj)
            % Get charges of the species in the chemical system
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     charges (float): Charges of the species in the chemical system

            charges = -obj.stoichiometricMatrix(:, obj.ind_E);
        end

        function [chargeIons, indexIons] = getChargeIons(obj)
            % Get charges of the ion species in the chemical system
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     chargeIons (float): Charges of the ion species in the chemical system

            % Default
            indexIons = [];

            % Get charges of the species
            charges = getCharges(obj);
            chargeIons = charges( isIonized(obj) );

            % Get index of ions
            if nargout > 1
                indexIons = obj.getIndexIons(obj.listSpecies);
            end

        end

    end

    methods (Static)

        function clearThermoCache()
            % Clear the cached species thermo curve handles used by evaluateSpeciesThermo
            combustiontoolbox.core.ChemicalSystem.speciesThermoCurveCache();
        end

    end

    methods (Access = private)
        
        function obj = setContainedElements(obj)
            % Obtain containted elements from the given set of species (reactants and products)
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the list of elements contained in the system
        
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

            % Get index of elements (equivalence ratio)
            obj.ind_C = find(ismember(obj.listElements, 'C'));
            obj.ind_H = find(ismember(obj.listElements, 'H'));
            obj.ind_O = find(ismember(obj.listElements, 'O'));
            obj.ind_N = find(ismember(obj.listElements, 'N'));
            obj.ind_E = find(ismember(obj.listElements, 'E'));
            obj.ind_S = find(ismember(obj.listElements, 'S'));
            obj.ind_Si = find(ismember(obj.listElements, 'SI'));
            obj.ind_B = find(ismember(obj.listElements, 'B'));
        end

        function index = getIndexIons(obj, species)
            % Get index of ions for the given list of chemical species
            %
            % Args:
            %     species (cell): List of chemical species
            %
            % Returns:
            %     index (float): Index of ions
        
            index = find( isIonized(obj, species) );
        end

        function obj = setIndexPhaseSpecies(obj)
            % Get index of gaseous, condensed and cryogenic species
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the index of gaseous, condensed and cryogenic species
        
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

            % Get index of species
            obj.indexSpecies = [obj.indexGas, obj.indexCondensed];
        end        

        function obj = sortIndexPhaseSpecies(obj)
            % Reorginize index of gaseous, condensed and cryogenic species
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %     LS (cell): Name list species / list of species
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the index of gaseous, condensed and cryogenic species sorted
            
            % Initialization
            obj.indexGas = []; obj.indexCondensed = []; obj.indexCryogenic = [];

            % Get index of gaseous, condensed and cryogenic species
            obj = obj.setIndexPhaseSpecies();
            
            % Get index of species
            obj.indexSpecies = [obj.indexGas, obj.indexCondensed];
        end

        function obj = setStoichiometricMatrix(obj)
            % Set stoichiometric matrix of the chemical system
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with the stoichiometric matrix set

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

        function obj = setStaticSpeciesProperties(obj)
            % Set static species properties from Species objects
            %
            % Args:
            %     obj (ChemicalSystem): ChemicalSystem object
            %
            % Returns:
            %     obj (ChemicalSystem): ChemicalSystem object with formationEnthalpy, formationInternalEnergy,
            %         molecularWeight, phase, temperatureMin, and temperatureMax set

            obj.formationEnthalpy = zeros(obj.numSpecies, 1);
            obj.formationInternalEnergy = zeros(obj.numSpecies, 1);
            obj.molecularWeight = zeros(obj.numSpecies, 1);
            obj.phase = zeros(obj.numSpecies, 1);
            obj.temperatureMin = zeros(obj.numSpecies, 1);
            obj.temperatureMax = zeros(obj.numSpecies, 1);

            for i = 1:obj.numSpecies
                species = obj.species.(obj.listSpecies{i});
                obj.formationEnthalpy(i) = species.hf;       % [J/mol]
                obj.formationInternalEnergy(i) = species.ef; % [J/mol]
                obj.molecularWeight(i) = species.W;          % [kg/mol]
                obj.phase(i) = species.phase;                % [-]
                obj.temperatureMin(i) = species.T(1);         % [K]
                obj.temperatureMax(i) = species.T(end);       % [K]
            end
        end

    end

    methods (Access = private, Static)

        function [cpCurves, h0Curves, s0Curves, g0Curves] = speciesThermoCurveCache(obj, index, mode)
            % Cache requested species thermo curve handles across stateless evaluations

            persistent cachedListSpecies cachedCPcurves cachedH0curves cachedS0curves cachedG0curves
            persistent FLAG_CACHE_CP FLAG_CACHE_H FLAG_CACHE_S FLAG_CACHE_G

            if nargin == 0
                cachedListSpecies = {};
                cachedCPcurves = {};
                cachedH0curves = {};
                cachedS0curves = {};
                cachedG0curves = {};
                FLAG_CACHE_CP = false;
                FLAG_CACHE_H = false;
                FLAG_CACHE_S = false;
                FLAG_CACHE_G = false;
                cpCurves = {};
                h0Curves = {};
                s0Curves = {};
                g0Curves = {};
                return
            end

            listSpecies = obj.listSpecies(:).';

            if isempty(cachedListSpecies) || ~isequal(cachedListSpecies, listSpecies)
                cachedListSpecies = listSpecies;
                numSystemSpecies = obj.numSpecies;
                cachedCPcurves = cell(1, numSystemSpecies);
                cachedH0curves = cell(1, numSystemSpecies);
                cachedS0curves = cell(1, numSystemSpecies);
                cachedG0curves = cell(1, numSystemSpecies);
                FLAG_CACHE_CP = false;
                FLAG_CACHE_H = false;
                FLAG_CACHE_S = false;
                FLAG_CACHE_G = false;
            end

            cpCurves = {};
            h0Curves = {};
            s0Curves = {};
            g0Curves = {};

            switch mode
                case 1 % Full thermo cache: cp, h, s, g
                    FLAG_UPDATE_CP = ~FLAG_CACHE_CP;
                    FLAG_UPDATE_H = ~FLAG_CACHE_H;
                    FLAG_UPDATE_S = ~FLAG_CACHE_S;
                    FLAG_UPDATE_G = ~FLAG_CACHE_G;

                    if FLAG_UPDATE_CP || FLAG_UPDATE_H || FLAG_UPDATE_S || FLAG_UPDATE_G
                        cacheMode = FLAG_UPDATE_CP + 2 * FLAG_UPDATE_H + 4 * FLAG_UPDATE_S + 8 * FLAG_UPDATE_G;

                        switch cacheMode
                            case 1 % Update cp
                                for i = 1:obj.numSpecies
                                    cachedCPcurves{i} = obj.species.(listSpecies{i}).cpcurve;
                                end
                            case 2 % Update h
                                for i = 1:obj.numSpecies
                                    cachedH0curves{i} = obj.species.(listSpecies{i}).h0curve;
                                end
                            case 3 % Update cp, h
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedH0curves{i} = species.h0curve;
                                end
                            case 4 % Update s
                                for i = 1:obj.numSpecies
                                    cachedS0curves{i} = obj.species.(listSpecies{i}).s0curve;
                                end
                            case 5 % Update cp, s
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedS0curves{i} = species.s0curve;
                                end
                            case 6 % Update h, s
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedH0curves{i} = species.h0curve;
                                    cachedS0curves{i} = species.s0curve;
                                end
                            case 7 % Update cp, h, s
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedH0curves{i} = species.h0curve;
                                    cachedS0curves{i} = species.s0curve;
                                end
                            case 8 % Update g
                                for i = 1:obj.numSpecies
                                    cachedG0curves{i} = obj.species.(listSpecies{i}).g0curve;
                                end
                            case 9 % Update cp, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                            case 10 % Update h, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedH0curves{i} = species.h0curve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                            case 11 % Update cp, h, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedH0curves{i} = species.h0curve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                            case 12 % Update s, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedS0curves{i} = species.s0curve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                            case 13 % Update cp, s, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedS0curves{i} = species.s0curve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                            case 14 % Update h, s, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedH0curves{i} = species.h0curve;
                                    cachedS0curves{i} = species.s0curve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                            case 15 % Update cp, h, s, g
                                for i = 1:obj.numSpecies
                                    species = obj.species.(listSpecies{i});
                                    cachedCPcurves{i} = species.cpcurve;
                                    cachedH0curves{i} = species.h0curve;
                                    cachedS0curves{i} = species.s0curve;
                                    cachedG0curves{i} = species.g0curve;
                                end
                        end

                        FLAG_CACHE_CP = true;
                        FLAG_CACHE_H = true;
                        FLAG_CACHE_S = true;
                        FLAG_CACHE_G = true;
                    end

                    cpCurves = cachedCPcurves(index);
                    h0Curves = cachedH0curves(index);
                    s0Curves = cachedS0curves(index);
                    g0Curves = cachedG0curves(index);
                case 2 % HCPS cache: h, cp, s
                    FLAG_UPDATE_CP = ~FLAG_CACHE_CP;
                    FLAG_UPDATE_H = ~FLAG_CACHE_H;
                    FLAG_UPDATE_S = ~FLAG_CACHE_S;

                    if FLAG_UPDATE_CP || FLAG_UPDATE_H || FLAG_UPDATE_S
                        if FLAG_UPDATE_CP && FLAG_UPDATE_H && FLAG_UPDATE_S
                            for i = 1:obj.numSpecies
                                species = obj.species.(listSpecies{i});
                                cachedCPcurves{i} = species.cpcurve;
                                cachedH0curves{i} = species.h0curve;
                                cachedS0curves{i} = species.s0curve;
                            end
                        elseif FLAG_UPDATE_CP && FLAG_UPDATE_H
                            for i = 1:obj.numSpecies
                                species = obj.species.(listSpecies{i});
                                cachedCPcurves{i} = species.cpcurve;
                                cachedH0curves{i} = species.h0curve;
                            end
                        elseif FLAG_UPDATE_CP && FLAG_UPDATE_S
                            for i = 1:obj.numSpecies
                                species = obj.species.(listSpecies{i});
                                cachedCPcurves{i} = species.cpcurve;
                                cachedS0curves{i} = species.s0curve;
                            end
                        elseif FLAG_UPDATE_H && FLAG_UPDATE_S
                            for i = 1:obj.numSpecies
                                species = obj.species.(listSpecies{i});
                                cachedH0curves{i} = species.h0curve;
                                cachedS0curves{i} = species.s0curve;
                            end
                        elseif FLAG_UPDATE_CP
                            for i = 1:obj.numSpecies
                                cachedCPcurves{i} = obj.species.(listSpecies{i}).cpcurve;
                            end
                        elseif FLAG_UPDATE_H
                            for i = 1:obj.numSpecies
                                cachedH0curves{i} = obj.species.(listSpecies{i}).h0curve;
                            end
                        else
                            for i = 1:obj.numSpecies
                                cachedS0curves{i} = obj.species.(listSpecies{i}).s0curve;
                            end
                        end

                        FLAG_CACHE_CP = true;
                        FLAG_CACHE_H = true;
                        FLAG_CACHE_S = true;
                    end

                    cpCurves = cachedCPcurves(index);
                    h0Curves = cachedH0curves(index);
                    s0Curves = cachedS0curves(index);
                case 3 % CPS cache: cp, s
                    FLAG_UPDATE_CP = ~FLAG_CACHE_CP;
                    FLAG_UPDATE_S = ~FLAG_CACHE_S;

                    if FLAG_UPDATE_CP || FLAG_UPDATE_S
                        if FLAG_UPDATE_CP && FLAG_UPDATE_S
                            for i = 1:obj.numSpecies
                                species = obj.species.(listSpecies{i});
                                cachedCPcurves{i} = species.cpcurve;
                                cachedS0curves{i} = species.s0curve;
                            end
                        elseif FLAG_UPDATE_CP
                            for i = 1:obj.numSpecies
                                cachedCPcurves{i} = obj.species.(listSpecies{i}).cpcurve;
                            end
                        else
                            for i = 1:obj.numSpecies
                                cachedS0curves{i} = obj.species.(listSpecies{i}).s0curve;
                            end
                        end

                        FLAG_CACHE_CP = true;
                        FLAG_CACHE_S = true;
                    end

                    cpCurves = cachedCPcurves(index);
                    s0Curves = cachedS0curves(index);
                case 4 % H cache: h
                    if ~FLAG_CACHE_H
                        for i = 1:obj.numSpecies
                            cachedH0curves{i} = obj.species.(listSpecies{i}).h0curve;
                        end

                        FLAG_CACHE_H = true;
                    end

                    h0Curves = cachedH0curves(index);
                case 5 % HG cache: h, g
                    FLAG_UPDATE_H = ~FLAG_CACHE_H;
                    FLAG_UPDATE_G = ~FLAG_CACHE_G;

                    if FLAG_UPDATE_H || FLAG_UPDATE_G
                        if FLAG_UPDATE_H && FLAG_UPDATE_G
                            for i = 1:obj.numSpecies
                                species = obj.species.(listSpecies{i});
                                cachedH0curves{i} = species.h0curve;
                                cachedG0curves{i} = species.g0curve;
                            end
                        elseif FLAG_UPDATE_H
                            for i = 1:obj.numSpecies
                                cachedH0curves{i} = obj.species.(listSpecies{i}).h0curve;
                            end
                        else
                            for i = 1:obj.numSpecies
                                cachedG0curves{i} = obj.species.(listSpecies{i}).g0curve;
                            end
                        end

                        FLAG_CACHE_H = true;
                        FLAG_CACHE_G = true;
                    end

                    h0Curves = cachedH0curves(index);
                    g0Curves = cachedG0curves(index);
                case 6 % G cache: g
                    if ~FLAG_CACHE_G
                        for i = 1:obj.numSpecies
                            cachedG0curves{i} = obj.species.(listSpecies{i}).g0curve;
                        end

                        FLAG_CACHE_G = true;
                    end

                    g0Curves = cachedG0curves(index);
                otherwise
                    error('Unknown species thermo cache mode.');
            end
        end

    end

end