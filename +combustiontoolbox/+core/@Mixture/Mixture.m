classdef Mixture < handle & matlab.mixin.Copyable
    % The :mat:class:`Mixture` class is used to store the properties of a chemical mixture.
    %
    % The :mat:class:`Mixture` object can be initialized as follows: ::
    %
    %      mix = Mixture(chemicalSystem)
    %
    % This creates an instance of the `Mixture` class and initializes it with a predefined chemical system.
    %
    % See also: :mat:func:`ChemicalSystem`, :mat:func:`Database`, :mat:func:`NasaDatabase`

    properties
        T (1,1)               % Temperature [K]
        p (1,1)               % Pressure [bar]
        N                     % Total number of moles [mol]
        hf                    % Enthalpy of formation [J]
        ef                    % Internal energy of formation [J]
        h                     % Enthalpy [J]
        e                     % Internal energy [J]
        g                     % Gibbs energy [J] 
        s                     % Entropy [J/K]
        cp                    % Specific heat at constant pressure [J/K]
        cv                    % Specific heat at constant volume [J/K]
        gamma                 % Adiabatic index [-]
        gamma_s               % Adiabatic index [-]
        sound                 % Speed of sound [m/s]
        s0                    % Entropy (frozen) [J/K]
        DhT                   % Thermal enthalpy [J]
        DeT                   % Thermal internal energy [J]
        Ds                    % Entropy of mixing [J/K]
        rho                   % Density [kg/m3]
        v                     % Volume [m3]
        vSpecific (1,1)       % Specific volume [m3/kg]
        W                     % Molecular weight [kg/mol]
        MW                    % Mean molecular weight [kg/mol]
        mi                    % Mass mixture [kg]
        Xi                    % Molar fractions [-]
        Yi                    % Mass fractions [-]
        phase                 % Phase vector [-]
        dVdT_p                % Derivative of volume with respect to temperature at constant pressure [-]
        dVdp_T                % Derivative of volume with respect to pressure at constant temperature [-]
        equivalenceRatio      % Equivalence ratio [-]
        equivalenceRatioSoot  % Theoretical equivalence ratio at which soot may appear [-]
        stoichiometricMoles   % Theoretical moles of the oxidizer of reference for a stoichiometric combustion
        percentageFuel        % Percentage of fuel in the mixture [%]
        fuelOxidizerMassRatio % Mass ratio of oxidizer to fuel [-]
        oxidizerFuelMassRatio % Mass ratio of fuel to oxidizer [-]
        natomElements         % Vector atoms of each element [-]
        natomElementsReact    % Vector atoms of each element without frozen species [-]
        chemicalSystem        % Chemical system object
        equationState         % Equation of State object
        u                     % Velocity module relative to the shock front [m/s]
        uShock                % Velocity module in the shock tube [m/s]
        uNormal               % Normal component of u [m/s]
        cjSpeed               % Chapman-Jouguet speed
        mach                  % Mach number [-]
        driveFactor           % Overdriven/Underdriven factor (detonations)
        beta                  % Wave angle [deg]
        theta                 % Deflection angle [deg]
        betaMin               % Minimum wave angle [deg]
        betaMax               % Maximum wave angle [deg]
        thetaMin              % Minimum eflection angle [deg]
        thetaMax              % Maximum deflection angle [deg]
        betaSonic             % Wave angle at the sonic point [deg]
        thetaSonic            % Deflection angle at the sonic point [deg]
        indexMin              % Index of the minimum deflection angle
        indexMax              % Index of the maximum deflection angle
        indexSonic            % Index of the sonic point
        polar                 % Properties of the polar solution
        areaRatio             % Area ratio = area_i / areaThroat
        areaRatioChamber      % Area ratio = areaChamber / areaThroat
        cstar                 % Characteristic velocity [m/s]
        cf                    % Thrust coefficient [-]
        I_sp                  % Specific impulse [s]
        I_vac                 % Vacuum impulse [s]
        config                % Mixture configuration object
    end

    properties (Access = private)
        indexSpecies           % Index of the species (initial mixture)
        indexGas               % Index of the gas species (initial mixture)
        Tspecies               % Species-specific initial temperatures [K] (initial mixture)
        FLAG_TSPECIES = false  % Flag to indicate species-specific initial temperatures are defined (initial mixture)
        FLAG_VOLUME = false    % Flag to indicate specific volume is defined (initial mixture)
    end
    
    properties (Hidden)
        errorMoles = 0         % Relative error in the moles calculation [-]
        errorMolesIons = 0     % Relative error in the moles of ions calculation [-]
        errorProblem = 0       % Relative error in the problem [-]
        cp_f                   % Frozen component of the specific heat at constant pressure
        cp_r                   % Reactive component of the specific heat at constant pressure
        dNi_T                  % Partial derivative of the number of moles with respect to temperature
        dN_T                   % Partial derivative of the total number of moles with respect to temperature
        dNi_p                  % Partial derivative of the number of moles with respect to pressure
        dN_p                   % Partial derivative of the total number of moles with respect to pressure
        chemicalSystemProducts % Chemical system containing only the list of products
        problemType            % Problem type
        rangeName              % Parametric property name
        quantity               % Composition (initial mixture)
        numSpecies             % Number of species (initial mixture)
        listSpecies            % List of species (initial mixture)
        listSpeciesFuel        % List of species fuel (initial mixture)
        listSpeciesOxidizer    % List of species oxidizer (initial mixture)
        listSpeciesInert       % List of species inert (initial mixture)
        molesFuel              % Moles of fuel (initial mixture)
        molesOxidizer          % Moles of oxidizer (initial mixture)
        molesInert             % Moles of inert (initial mixture)
        ratioOxidizer          % Ratio oxidizer relative to the oxidizer of reference (initial mixture)
        fuel                   % Fuel atoms (initial mixture)
        FLAG_FUEL = false      % Flag to indicate fuel species are defined (initial mixture)
        FLAG_OXIDIZER = false  % Flag to indicate oxidizer species are defined (initial mixture)
        FLAG_INERT = false     % Flag to indicate inert species are defined (initial mixture)
        FLAG_REACTION = false  % Flag to indicate chemical reaction is defined
    end
    
    methods
        
        mix = setStagnation(mix, varargin)
        print(mix, varargin)

        function obj = Mixture(chemicalSystem, varargin)
            % Mixture constructor

            % Definitions
            defaultTemperature = 300; % [K]
            defaultPressure = 1;      % [bar]
            defaultEoS = combustiontoolbox.core.EquationStateIdealGas();
            defaultConfig = combustiontoolbox.core.MixtureConfig();

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'chemicalSystem', @(x) isa(x, 'combustiontoolbox.core.ChemicalSystem'));
            addOptional(ip, 'T', defaultTemperature, @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'p', defaultPressure, @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'eos', defaultEoS, @(x) isa(x, 'combustiontoolbox.core.EquationState'));
            addParameter(ip, 'config', defaultConfig, @(x) isa(x, 'combustiontoolbox.core.MixtureConfig'));
            parse(ip, chemicalSystem, varargin{:});
            
            % Set properties
            obj.chemicalSystem = ip.Results.chemicalSystem;
            obj.T = ip.Results.T;
            obj.p = ip.Results.p;
            obj.equationState = ip.Results.eos;
            obj.config = ip.Results.config;
        end

        function obj = setTemperature(obj, T, varargin)
            % Set temperature [K] and compute thermodynamic properties
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     T (float): Temperature [K]
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            %
            % Example:
            %     setTemperature(obj, 300)

            % Definitions
            defaultUnits = 'K';

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'T', @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'units', defaultUnits, @(x) ischar(x));
            parse(ip, T, varargin{:});
            
            % Change units to Kelvin
            if ~strcmpi(ip.Results.units, 'K')
                T = combustiontoolbox.common.Units.convert(T, ip.Results.units, 'K');
            end

            % Assign temperature
            obj.T = T;
            
            % Update thermodynamic state
            updateThermodynamics(obj);
        end

        function obj = setPressure(obj, p, varargin)
            % Set pressure [bar] and compute thermodynamic properties
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     p (float): Pressure [bar]
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            %
            % Example:
            %     setPressure(obj, 1)
            
            % Definitions
            defaultUnits = 'bar';

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'p', @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'units', defaultUnits, @(x) ischar(x));
            parse(ip, p, varargin{:});
            
            % Change units to bar
            if ~strcmpi(ip.Results.units, 'bar')
                p = Units.convert(p, ip.Results.units, 'bar');
            end

            % Assign pressure
            obj.p = p;

            % Update thermodynamic state
            updateThermodynamics(obj);
        end

        function obj = setVolume(obj, vSpecific, varargin)
            % Set specific volume [m3/kg] and compute thermodynamic properties
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     vSpecific (float): Specific volume [m3/kg]
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            %
            % Example:
            %     setVolume(obj, 1)

            % Definitions
            defaultUnits = 'm3/kg';

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'vSpecific', @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'units', defaultUnits, @(x) ischar(x));
            parse(ip, vSpecific, varargin{:});
            
            % Change units to m3/kg
            if ~strcmpi(ip.Results.units, 'm3/kg')
                vSpecific = combustiontoolbox.common.Units.convert(vSpecific, ip.Results.units, 'm3/kg');
            end

            % Assign specific volume
            obj.vSpecific = vSpecific;
            obj.FLAG_VOLUME = true;
            
            % Update thermodynamic state
            updateThermodynamics(obj);
        end

        function obj = set(obj, listSpecies, varargin)
            % Set species and quantity and compute thermodynamic properties
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     listSpecies (cell): List of species
            %
            % Optional Args:
            %     type (char): Type of species (fuel, oxidizer, inert)
            %     quantity (float): Quantity of species [mol, weightPercentage]
            %     units (char): Units of quantity (mol, weightPercentage)
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            %
            % Examples:
            %     * set(obj, {'H2'}, 1)
            %     * set(obj, {'H2'}, 'fuel', 1)
            %     * set(obj, {'CH6N2bLb', 'N2H4bLb'}, 'fuel', [86, 14], 'weightPercentage');
            %     * set(obj, {'N2', 'O2'}, 'oxidizer', [79, 21]);

            % Import packages
            import combustiontoolbox.utils.findIndex
            
            % Definitions
            quantityDefault = 1;
            unitsDefault = 'mol'; % mol or weightPercentage

            % 
            if nargin > 3
                type = varargin{1};
                quantity = varargin{2};

                if nargin > 4
                    unitsDefault = varargin{3};
                    
                    switch lower(unitsDefault)
                        case 'mol'
                            % Nothing to do
                        case 'weightpercentage'
                            quantity = combustiontoolbox.common.Units.convertWeightPercentage2moles(listSpecies, quantity, obj.chemicalSystem.database);
                    end

                end

                switch lower(type)
                    case 'fuel'
                        obj.listSpeciesFuel = [obj.listSpeciesFuel, listSpecies];
                        obj.molesFuel = [obj.molesFuel, quantity];
                        obj.FLAG_FUEL = true;
                    case 'oxidizer'
                        obj.listSpeciesOxidizer = [obj.listSpeciesOxidizer, listSpecies];
                        obj.molesOxidizer = [obj.molesOxidizer, quantity];
                        if isempty(obj.ratioOxidizer), obj.ratioOxidizer = obj.molesOxidizer; end
                        obj.FLAG_OXIDIZER = true;
                    case 'inert'
                        obj.listSpeciesInert = [obj.listSpeciesInert, listSpecies];
                        obj.molesInert = [obj.molesInert, quantity];
                        obj.FLAG_INERT = true;
                end

            elseif nargin > 2
                quantity = varargin{1};
            end

            % Update listSpecies, quantity and numSpecies of the initial mixture
            obj.listSpecies = [obj.listSpecies, listSpecies];
            obj.quantity = [obj.quantity, quantity];
            obj.numSpecies = length(obj.listSpecies);

            % Check if species are contained in the chemical system
            obj.chemicalSystem.checkSpecies(listSpecies);
            
            % Update index species in the mixture
            updateIndexSpecies(obj);

            % Get index react species
            obj.chemicalSystem.setReactIndex(obj.listSpeciesInert);
            
            % Get index products species
            obj.chemicalSystem.indexProducts = findIndex(obj.chemicalSystem.listSpecies, obj.chemicalSystem.listProducts);
            
            % Get index gas species
            obj.indexGas = find(ismember(obj.listSpecies, obj.chemicalSystem.listSpecies(obj.chemicalSystem.indexGas)));

            % Update composition
            updateComposition(obj);

            % Update thermodynamic state
            updateThermodynamics(obj);
        end

        function obj = setEquivalenceRatio(obj, equivalenceRatio)
            % Set equivalence ratio and compute thermodynamic properties
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     equivalenceRatio (float): Equivalence ratio [-]
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            %
            % Example:
            %     setEquivalenceRatio(obj, 1)

            % Definitions
            obj.equivalenceRatio = equivalenceRatio;
            
            % Update composition
            updateComposition(obj);

            % Compute thermodynamic properties
            updateThermodynamics(obj);
        end

        function obj = setTemperatureSpecies(obj, speciesTemperatures)
            % Set species-specific temperatures and update equilibrium temperature and properties
            %
            % Args:
            %     speciesTemperatures (float): Array with the temperatures [K] for each species in the mixture
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated equilibrium temperature and properties
            %
            % Note: The temperature vector is assigned to the species in the same order as they were added to the mixture object
            %
            % Example:
            %     setTemperatureSpecies(obj, [300, 400, 350])
            
            % Validate input: check that the number of temperatures equals the number of species
            if numel(speciesTemperatures) ~= obj.numSpecies
                error('Temperature input must be either a scalar or a vector of length equal to the number of species (%d).', obj.numSpecies);
            end       
            
            % Definitions
            obj.FLAG_TSPECIES = true;
            obj.Tspecies = speciesTemperatures;

            % Compute the equilibrium temperature
            obj.T = obj.computeEquilibriumTemperature(speciesTemperatures);
        end

        function obj = computeEquivalenceRatio(obj)
            % Compute equivalence ratio [-]
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated equivalence ratio [-]
            
            % Check if initial state is defined (temperature, pressure, and composition)
            if ~obj.FLAG_FUEL || ~obj.FLAG_OXIDIZER
                return
            end

            % Set oxidizer of reference
            obj.chemicalSystem.setOxidizerReference(obj.listSpeciesOxidizer);
            
            % Computation of theoretical stoichiometricMoles
            obj.defineF();
            
            % Define oxidizer propertiesMatrix
            obj.defineO();
            
            % Compute percentage Fuel, Oxidizer/Fuel ratio and equivalence ratio
            obj.computeRatiosFuelOxidizer(obj.chemicalSystem.propertiesMatrixFuel, obj.chemicalSystem.propertiesMatrixOxidizer);
        end

        function obj = computeRatiosFuelOxidizer(obj, propertiesMatrixFuel, propertiesMatrixOxidizer)
            % Compute percentage Fuel, Oxidizer/Fuel ratio and equivalence ratio
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     propertiesMatrixFuel (float): Properties matrix of the fuel
            %     propertiesMatrixOxidizer (float): Properties matrix of the oxidizer
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties

            if obj.FLAG_FUEL && obj.FLAG_OXIDIZER
                mass_fuel = obj.getMass(obj.chemicalSystem, propertiesMatrixFuel);
                mass_oxidizer = obj.getMass(obj.chemicalSystem, propertiesMatrixOxidizer);
                mass_mixture = obj.mi;
                obj.percentageFuel = mass_fuel / mass_mixture * 100;
                obj.oxidizerFuelMassRatio = mass_oxidizer / mass_fuel;
                obj.fuelOxidizerMassRatio = 1 / obj.oxidizerFuelMassRatio;
                FO_moles = sum(propertiesMatrixFuel(:, obj.chemicalSystem.ind_ni)) / sum(propertiesMatrixOxidizer(obj.chemicalSystem.oxidizerReferenceIndex, obj.chemicalSystem.ind_ni));
                FO_moles_st = abs(sum(propertiesMatrixFuel(:, obj.chemicalSystem.ind_ni)) / (obj.fuel.C + obj.fuel.H / 4 - obj.fuel.O / 2 + obj.fuel.S + obj.fuel.Si + 3/4 * obj.fuel.B) * (0.5 * obj.chemicalSystem.oxidizerReferenceAtomsO));
                obj.equivalenceRatio = FO_moles / FO_moles_st;
                computeEquivalenceRatioSoot(obj); 
            elseif obj.FLAG_FUEL
                obj.percentageFuel = 100;
                obj.fuelOxidizerMassRatio = inf;
                obj.oxidizerFuelMassRatio = 0;
                obj.equivalenceRatio = [];
                obj.equivalenceRatioSoot = [];
            else
                obj.percentageFuel = 0;
                obj.fuelOxidizerMassRatio = 0;
                obj.oxidizerFuelMassRatio = inf;
                obj.equivalenceRatio = [];
                obj.equivalenceRatioSoot = [];
            end
        
        end

        function obj = computeEquivalenceRatioSoot(obj)
            % Compute guess of equivalence ratio in which soot appears considering complete combustion
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with theoretical equivalence ratio at which soot appears [-]
        
            obj.equivalenceRatioSoot = 2 / (obj.fuel.C - obj.fuel.O) * (obj.fuel.C + obj.fuel.H / 4 - obj.fuel.O / 2);
        
            if obj.equivalenceRatioSoot <= 1e-5 || isnan(obj.equivalenceRatioSoot)
                obj.equivalenceRatioSoot = inf;
            end
        
        end

        function objArray = setProperties(obj, property, value, varargin)
            % Obtain properties at equilibrium for the given thermochemical transformation
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     property (char): Property to be set
            %     value (float): Value of the property
            %
            % Optional Args (key-value pairs):
            %     * property (char): Property to be set
            %     * value (float): Value of the property
            %
            % Returns:
            %     objArray (Mixture): Array of Mixture objects with the computed properties
            %
            % Note:
            %     * Use this method after setting the initial composition of the mixture
            %
            % Examples:
            %     * mixArray = setProperties(mix, 'equivalenceRatio', value)
            %     * mixArray = setProperties(mix, 'equivalenceRatio', value, 'temperature', value)
            %     * mixArray = setProperties(mix, 'equivalenceRatio', value, 'temperature', value, 'pressure', value)
            %     * mixArray = setProperties(mix, 'phi', value, 'T', value, 'p', value)
            
            % Initialization
            FLAG_MACH = false;
            
            % Assign value to the property
            properties = {property, varargin{1:2:end}};
            values = {value, varargin{2:2:end}};
            
            % Definitions
            numProperties = min(length(properties), length(values));

            % Check vectors
            FLAG_VECTOR = cellfun(@(x) numel(x) > 1, values);
            
            if sum(FLAG_VECTOR)
                % Create vectors same size
                FLAG_VECTOR_FIRST = find(FLAG_VECTOR, 1);
                aux = ones(size(values{FLAG_VECTOR_FIRST}));
    
                for i = find(~FLAG_VECTOR)
                    values(i) = {values{i} * aux};
                end
                
                % Define property parametric study
                obj.rangeName = properties{FLAG_VECTOR_FIRST};
                
                % Check rangeName to match Mixture's property
                if any(~strcmpi({'T', 'p', 'vspecific', 'phi', 'u', 'mach', 'beta', 'theta', 'drive_factor', 'aratio', 'aratio_c'}, obj.rangeName))

                    switch lower(obj.rangeName)
                        case {'temperature'}
                            obj.rangeName = 'T';
                        case {'pressure'}
                            obj.rangeName = 'p';
                        case {'volume', 'v'}
                            obj.rangeName = 'vSpecific';
                        case {'phi'}
                            obj.rangeName = 'equivalenceRatio';
                        case {'velocity', 'u1'}
                            obj.rangeName = 'u';
                        case {'m1'}
                            obj.rangeName = 'mach';
                        case {'wave angle', 'waveangle', 'wave'}
                            obj.rangeName = 'beta';
                        case {'deflection angle', 'deflectionangle', 'deflection'}
                            obj.rangeName = 'theta';
                        case {'drive_factor'}
                            obj.rangeName = 'driveFactor';
                        case {'aratio'}
                            obj.rangeName = 'areaRatio';
                        case {'aratio_c'}
                            obj.rangeName = 'areaRatioChamber';
                    end

                end

            else
                FLAG_VECTOR_FIRST = 1;
            end

            % Get number of cases
            numCases = length(values{FLAG_VECTOR_FIRST});
            
            % Initialization
            objArray = obj.empty(0, numCases);

            % Set properties
            for j = 1:numCases
                % Create a copy of the mixture
                objArray(j) = obj.copyDeep();

                % Set property
                for i = 1:numProperties

                    switch lower(properties{i})
                        case {'temperature', 't'}
                            objArray(j).T = values{i}(j);
                        case {'pressure', 'p'}
                            objArray(j).p = values{i}(j);
                            objArray(j).FLAG_VOLUME = false;
                        case {'volume', 'vspecific', 'v'}
                            objArray(j).vSpecific = values{i}(j);
                            objArray(j).FLAG_VOLUME = true;
                        case {'equivalenceratio', 'phi'}
                            objArray(j).equivalenceRatio = values{i}(j);
                            objArray(j).updateComposition();
                        case {'velocity', 'u', 'u1'}
                            objArray(j).u = values{i}(j);
                        case {'mach', 'm1'}
                            objArray(j).mach = values{i}(j);
                            FLAG_MACH = true;
                        case {'wave angle', 'waveangle', 'wave', 'beta'}
                            objArray(j).beta = values{i}(j);
                        case {'deflection angle', 'deflectionangle', 'deflection', 'theta'}
                            objArray(j).theta = values{i}(j);
                        case {'drive_factor', 'drivefactor'}
                            objArray(j).driveFactor = values{i}(j);
                        case {'arearatio', 'aratio'}
                            objArray(j).areaRatio = values{i}(j);
                        case {'arearatiochamber', 'aratio_c'}
                            objArray(j).areaRatioChamber = values{i}(j);
                        otherwise
                            error('Property not found');
                    end

                end

                % Compute thermodynamic state of the mixture
                objArray(j).updateThermodynamics();

                % Additional inputs
                if FLAG_MACH
                    objArray(j).u = objArray(j).mach * objArray(j).sound;
                end

            end

        end

        function obj = updateThermodynamics(obj)
            % Update the thermodynamic state of the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            
            % Check if initial state is defined (temperature, pressure/volume, and composition)
            if ~sum(obj.quantity) && ~obj.T || (~obj.p && ~obj.vSpecific)
                return
            end

            if obj.FLAG_VOLUME
                % Compute molar volume [m3/mol] from specific volume [m3/kg]
                vMolar = vSpecific2vMolar(obj, obj.vSpecific, obj.quantity, obj.quantity(obj.indexGas));
                % Compute pressure in Pascals using the equationState
                pressure = obj.equationState.getPressure(obj.T, vMolar, obj.chemicalSystem.listSpecies, obj.quantity / sum(obj.quantity)); % [Pa]
                % Convert pressure to [bar]
                obj.p = pressure * combustiontoolbox.common.Units.Pa2bar;
            end
            
            % Assign values to the propertiesMatrix
            obj.chemicalSystem.setPropertiesMatrixInitialIndex(obj.listSpecies, obj.quantity, obj.T, obj.indexSpecies);

            % Compute thermodynamic properties
            computeThermodynamics(obj);
        end

        function obj = updateComposition(obj)
            % Update the composition of the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated properties
            
            % Check if initial composition is defined
            if ~sum(obj.quantity)
                return
            end

            % Check if mixture is compound of a fuel and an oxidizer
            if obj.FLAG_FUEL && obj.FLAG_OXIDIZER
                % Set oxidizer of reference
                obj.chemicalSystem.setOxidizerReference(obj.listSpeciesOxidizer);
                
                % Computation of theoretical stoichiometricMoles
                obj.defineF();
                
                % Define moles Oxidizer
                if ~isempty(obj.equivalenceRatio)
                    if isempty(obj.ratioOxidizer), obj.ratioOxidizer = obj.molesOxidizer; end
                    obj.molesOxidizer = obj.stoichiometricMoles / obj.equivalenceRatio .* obj.ratioOxidizer;
                end

                % Define oxidizer propertiesMatrix
                obj.defineO();
    
                % Update quantity
                obj.quantity = [obj.molesFuel, obj.molesOxidizer, obj.molesInert];
            end
            
            % Assign values to the propertiesMatrix
            obj.chemicalSystem.setPropertiesMatrixCompositionInitialIndex(obj.listSpecies, obj.quantity, obj.indexSpecies);

            % Compute composition
            computeComposition(obj);

            % Compute equivalence ratio
            computeEquivalenceRatio(obj);

            % Check complete combustion
            if ~isempty(obj.equivalenceRatio)
                checkCompleteReaction(obj.chemicalSystem, obj.equivalenceRatio, obj.equivalenceRatioSoot);
            end
            
            % Get system containing only the list of products
            obj.chemicalSystemProducts = getSystemProducts(obj.chemicalSystem);
        end

        function obj = updateIndexSpecies(obj)
            % Update index species in the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated indexSpecies

            obj.indexSpecies = combustiontoolbox.utils.findIndex(obj.chemicalSystem.listSpecies, obj.listSpecies);
        end

        function vMolar = vSpecific2vMolar(obj, vSpecific, moles, molesGas, varargin)
            % Compute molar volume [m3/mol] from specific volume [m3/kg]

            % Get index specie
            if nargin == 4
                index = combustiontoolbox.utils.findIndex(obj.chemicalSystem.listSpecies, obj.listSpecies);
            else
                index = varargin{1};
            end
            
            % Compute Mean Molecular Weight [kg/mol]
            MW = computeMeanMolecularWeight(obj, moles, index);

            % Compute specific volume [m3/mol]
            vMolar = vSpecific * MW * sum(moles) / sum(molesGas);
        end

        function typeSpecies = getTypeSpecies(obj)
            % Create cell array with the type of species in the mixture
            %
            % Args:
            %     obj (Mixture): Mixture class
            %
            % Returns:
            %     typeSpecies (cell): Cell array with the type of species in the mixture
        
            typeFuel = repmat({'Fuel'}, size(obj.listSpeciesFuel));
            typeOxidizer = repmat({'Oxidizer'}, size(obj.listSpeciesOxidizer));
            typeInert = repmat({'Inert'}, size(obj.listSpeciesInert));
            typeSpecies = [typeInert, typeOxidizer, typeFuel];
        end
        
    end
    
    methods(Access = protected)
      
      function objCopy = copyElement(obj)
         % Override copyElement method:

         % Make a shallow copy of all properties
         objCopy = copyElement@matlab.mixin.Copyable(obj);
      end

      function objCopy = copyDeep(obj)
         % Make a deep copy of obj

         % Make a shallow copy of all properties
         objCopy = copyElement(obj);

         % Make a deep copy of the ChemicalSystem object
         objCopy.chemicalSystem = obj.chemicalSystem.copy();
      end
      
    end

    methods (Access = private)

        % function obj = plus(obj, obj2)
        %     % Merge two mixtures (operator overloading)

        %     % Define a common chemical system
            
        %     % Merge chemical composition

        %     % Compute common temperature and pressure
            
        %     % Compute properties of the mixture

        % end

        function obj = computeProperties(obj)
            % Compute composition and thermodynamic properties of the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with the computed properties
            %
            % Example:
            %     mix = computeProperties(obj)

            % Compute composition
            computeComposition(obj);
            
            % Compute thermodynamic properties
            computeThermodynamics(obj);
        end

        function computeComposition(obj)
            % Compute the composition of the mixture

            % Definitions
            system = obj.chemicalSystem;
            propertiesMatrix = system.propertiesMatrix; % Properties matrix

            % Unpack propertiesMatrix
            Ni = propertiesMatrix(:, system.ind_ni); % [mol]
            obj.N = sum(propertiesMatrix(:, system.ind_ni)); % [mol]
            obj.phase = propertiesMatrix(:, system.ind_phase); % [bool]

            % Compute total composition of gas species [mol]
            N_gas = sum(Ni(~obj.phase));

            % Compute molar fractions [-]
            obj.Xi = Ni / obj.N;

            % Compute molecular weight (gases) [kg/mol]
            obj.W = dot(Ni, propertiesMatrix(:, system.ind_W)) / N_gas;

            % Compute mean molecular weight [kg/mol]
            obj.MW = dot(Ni, propertiesMatrix(:, system.ind_W)) / obj.N;

            % Compute mass mixture [kg]
            obj.mi = obj.MW * obj.N; % [kg]

            % Compute mass fractions [-]
            obj.Yi = (Ni .* propertiesMatrix(:, system.ind_W)) ./ obj.mi;

            % Compute vector atoms of each element
            obj.natomElements = sum(Ni .* system.stoichiometricMatrix, 1);

            % Compute vector atoms of each element without frozen species
            obj.natomElementsReact = sum(propertiesMatrix(system.indexReact, system.ind_ni) .* system.stoichiometricMatrix(system.indexReact, :), 1);
        end

        function obj = computeThermodynamics(obj)
            % Compute thermodynamic properties of the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object with the computed properties
            %
            % Example:
            %     mix = computeThermodynamics(obj)

            if obj.FLAG_TSPECIES
                obj.setTemperatureSpecies(obj.Tspecies);
                % Reset species-specific temperatures 
                obj.Tspecies = [];
                obj.FLAG_TSPECIES = false;
            end
            
            % Definitions
            temperature = obj.T;
            pressure = obj.p;
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
            system = obj.chemicalSystem;
            propertiesMatrix = system.propertiesMatrix; % Properties matrix
            
            % Unpack propertiesMatrix
            Ni = propertiesMatrix(:, system.ind_ni); % [mol]
            obj.hf = dot(propertiesMatrix(:, system.ind_hfi), Ni); % [J]
            obj.h = dot(propertiesMatrix(:, system.ind_hi), Ni); % [J]
            obj.ef = dot(propertiesMatrix(:, system.ind_efi), Ni); % [J]
            obj.cp = dot(propertiesMatrix(:, system.ind_cpi), Ni); % [J/K]
            obj.s0 = dot(propertiesMatrix(:, system.ind_si), Ni); % [J/K]

            % Compute total composition of gas species [mol]
            N_gas = sum(Ni(~obj.phase));

            % Get non zero species
            FLAG_NONZERO = obj.Xi > 0;
            
            % Compute volume [m3]
            if N_gas
                obj.v = obj.equationState.getVolume(temperature, pressure * combustiontoolbox.common.Units.bar2Pa, obj.chemicalSystem.listSpecies, obj.Xi) * N_gas;
            else
                % Mixture that only has condensed species
                obj.v = obj.vSpecific * obj.mi;
            end

            % Compute specific volume [m3/kg]
            obj.vSpecific = obj.v / obj.mi;

            % Compute density [kg/m3]
            obj.rho = 1 / obj.vSpecific;

            % Compute internal energy [J]
            obj.e = obj.h - N_gas * R0 * temperature;

            % Compute thermal internal energy [J]
            obj.DeT = obj.e - obj.ef;

            % Compute thermal enthalpy [J]
            obj.DhT = obj.h - obj.hf;

            % Compute entropy of mixing [J/K]
            obj.Ds = obj.computeEntropyMixing(Ni, N_gas, R0, FLAG_NONZERO);

            % Compute entropy [J/K]
            obj.s = obj.s0 + obj.Ds;
            
            % Compute Gibbs energy [J]
            obj.g = obj.h - obj.T * obj.s;
            
            % Compute thermodynamic derivatives, cp, cv, gamma, and speed
            % of sound considering chemical reaction
            if obj.FLAG_REACTION
                % Set thermodynamic derivatives
                obj.dVdT_p = obj.dN_T + 1; % [-]
                obj.dVdp_T = obj.dN_p - 1; % [-]
        
                if ~any(isnan(obj.dNi_T)) && ~any(isinf(obj.dNi_T))
                    % Definitions
                    delta = ~obj.phase;
                    h0_j = propertiesMatrix(:, system.ind_hi); % [J/mol]

                    % Compute specific heat at constant pressure [J/K]
                    obj.cp_r = sum(h0_j / temperature .* (1 + delta .* (Ni - 1)) .* obj.dNi_T, 'omitnan');
                    obj.cp_f = obj.cp;
                    obj.cp = obj.cp_f + obj.cp_r;

                    % Compute specific heat at constant volume [J/K]
                    obj.cv = obj.cp + (N_gas * R0 * obj.dVdT_p^2) / obj.dVdp_T;

                    % Compute Adibatic index [-]
                    obj.gamma = obj.cp / obj.cv;
                    obj.gamma_s = -obj.gamma / obj.dVdp_T;

                    % Compute sound velocity [m/s]
                    obj.sound = sqrt(obj.gamma_s * combustiontoolbox.common.Units.bar2Pa * pressure / obj.rho);

                    % Compute Mach number
                    if ~isempty(obj.u)
                        obj.mach = obj.u / obj.sound;
                    end

                    return
                end

                obj.gamma_s = -1 / obj.dVdp_T; % [-]
                
                return
            end
            
            % Compute thermodynamic derivatives, cp, cv, gamma, and speed
            % of sound considering frozen chemistry
            
            % Set thermodynamic derivatives [-]
            obj.dVdT_p = 1;
            obj.dVdp_T = -1;

            % Compute specific heat at constant volume [J/K]
            obj.cv = obj.cp - R0 * N_gas;

            % Compute Adibatic index [-]
            obj.gamma = obj.cp / obj.cv;
            obj.gamma_s = obj.gamma;

            % Compute sound velocity [m/s]
            obj.sound = sqrt(obj.gamma * combustiontoolbox.common.Units.bar2Pa * pressure / obj.rho);
            
            % Compute Mach number
            if ~isempty(obj.u)
                obj.mach = obj.u / obj.sound;
            end

        end

        function MW = computeMeanMolecularWeight(obj, moles, index)
            % Compute Mean Molecular Weight [kg/mol]
            MW = dot(moles, obj.chemicalSystem.propertiesMatrix(index, obj.chemicalSystem.ind_W)) / sum(moles);
        end

        function Ds = computeEntropyMixing(obj, Ni, N_gas, R0, FLAG_NONZERO)
            % Compute entropy of mixing [J/K]
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     Ni (float): Vector of moles of each species [mol]
            %     N_gas (float): Total moles of gas species [mol]
            %     R0 (float): Universal gas constant [J/(mol-K)]
            %     FLAG_NONZERO (bool): Vector of nonzero species
            %
            % Returns:
            %     DS (float): Entropy of mixing [J/K]
            %
            % Note: 
            %     only nonzero for gaseous species
            %
            % Example:
            %     Ds = computeEntropyMixing(mix, Ni, N_gas, R0, FLAG_NONZERO)
        
            Dsi = Ni(FLAG_NONZERO) .* log(Ni(FLAG_NONZERO) / N_gas * obj.p) .* (1 - obj.phase(FLAG_NONZERO));
            Ds = -R0 * sum(Dsi);
        end

        function obj = defineF(obj)
            % Set Fuel of the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object

            if obj.FLAG_FUEL
                % Create a deep copy for the fuel mixture
                mixFuel = obj.copyDeep();
                % Set temperature-dependent matrix properties to zero
                mixFuel.chemicalSystem.clean();
                % Fill properties matrix with only fuel species
                mixFuel.chemicalSystem.setPropertiesMatrixComposition(obj.listSpeciesFuel, obj.molesFuel);
                % Compute composition properties 
                mixFuel.computeComposition();
                % Assign values elements C, H, O, N, S, and Si
                obj.assignAtomElementsFuel(mixFuel.natomElements);
                % Compute theoretical moles of the oxidizer of reference for a stoichiometric combustion
                obj.stoichiometricMoles = abs(obj.fuel.C + obj.fuel.H / 4 - obj.fuel.O / 2 ...
                    + obj.fuel.S + obj.fuel.Si + 3/4 * obj.fuel.B) / (0.5 * mixFuel.chemicalSystem.oxidizerReferenceAtomsO);
                % Assign propertiesMatrixFuel
                obj.chemicalSystem.propertiesMatrixFuel = mixFuel.chemicalSystem.propertiesMatrix;
                return
            end

            obj.fuel.C = 0;
            obj.fuel.H = 0;
            obj.fuel.O = 0;
            obj.fuel.N = 0;
            obj.fuel.S = 0;
            obj.fuel.Si = 0;
            obj.fuel.B = 0;
            obj.stoichiometricMoles = 0;
        end

        function obj = defineO(obj)
            % Set Oxidizer of the mixture
            %
            % Args:
            %     obj (Mixture): Mixture object
            %
            % Returns:
            %     obj (Mixture): Mixture object

            if isempty(obj.listSpeciesOxidizer)
                return
            end

            % Create a copy of chemicalSystem
            system = obj.chemicalSystem.copy();
            % Set temperature-dependent matrix properties to zero
            system.clean();
            % Fill properties matrix with only oxidizer species
            system.setPropertiesMatrixComposition(obj.listSpeciesOxidizer, obj.molesOxidizer);
            % Assign propertiesMatrixOxidizer
            obj.chemicalSystem.propertiesMatrixOxidizer = system.propertiesMatrix;
        end

        function obj = assignAtomElementsFuel(obj, natomElementsFuel)
            % Assign atomic counts of the fuel
            %
            % Args:
            %     obj (Mixture): Mixture object
            %     natomElementsFuel(float): vector containing atomic counts for the fuel
            %
            % Returns:
            %     obj (Mixture): Mixture object with updated atomic counts for the fuel
            
            if isempty(obj.chemicalSystem.ind_C), obj.fuel.C = 0; else, obj.fuel.C = natomElementsFuel(obj.chemicalSystem.ind_C); end
            if isempty(obj.chemicalSystem.ind_H), obj.fuel.H = 0; else, obj.fuel.H = natomElementsFuel(obj.chemicalSystem.ind_H); end
            if isempty(obj.chemicalSystem.ind_O), obj.fuel.O = 0; else, obj.fuel.O = natomElementsFuel(obj.chemicalSystem.ind_O); end
            if isempty(obj.chemicalSystem.ind_N), obj.fuel.N = 0; else, obj.fuel.N = natomElementsFuel(obj.chemicalSystem.ind_N); end
            if isempty(obj.chemicalSystem.ind_S), obj.fuel.S = 0; else, obj.fuel.S = natomElementsFuel(obj.chemicalSystem.ind_S); end
            if isempty(obj.chemicalSystem.ind_Si), obj.fuel.Si = 0; else, obj.fuel.Si = natomElementsFuel(obj.chemicalSystem.ind_Si); end
            if isempty(obj.chemicalSystem.ind_B), obj.fuel.B = 0; else, obj.fuel.B = natomElementsFuel(obj.chemicalSystem.ind_B); end

        end

        function T = computeEquilibriumTemperature(obj, speciesTemperatures)
            % Compute the equilibrium temperature [K] for a mixture of n species, each initially at a specified
            % temperature. A Newton-Raphson method is used to iteratively adjust the temperature until equilibrium is reached
            %
            % Args:
            %     speciesTemperatures (float): Array containing the initial temperatures [K] for each species
            %
            % Returns:
            %     T (float): Equilibrium temperature [K] of the mixture
            %
            % Example:
            %     T = obj.computeEquilibriumTemperature([300, 400, 350])
            
            % Import packages
            import combustiontoolbox.core.Species.*

            % Convergence criteria and maximum number of iterations
            tol0 = 1e-3;
            itMax = 100;
            
            % If all species have the same temperature, return it directly
            if isscalar(unique(speciesTemperatures))
                T = speciesTemperatures(1);
                return;
            end

            % Check if condensed species can be only evaluated a particular temperature
            FLAG_FIXED = checkTemperatureSpecies(obj);
            if FLAG_FIXED
                T = max(speciesTemperatures);
                return
            end

            % Determine problem type (default to 'TP' if not set)
            if isempty(obj.problemType)
                problemType = 'TP';
            else
                problemType = obj.problemType;
            end
            
            % Choose the appropriate property functions based on problem type
            switch lower(problemType)
                case {'tv', 'ev', 'sv', 'volume', 'v'}
                    funCpOrCv = @getHeatCapacityVolume;
                    funHorE   = @getInternalEnergy;
                otherwise
                    funCpOrCv = @getHeatCapacityPressure;
                    funHorE   = @getEnthalpy;
            end
            
            % Evaluate properties at each species' initial temperature
            CpOrCv_0 = obj.getPropertyListSpecies(funCpOrCv, speciesTemperatures);
            HorE_0   = obj.getPropertyListSpecies(funHorE, speciesTemperatures);
            
            % Initial guess for the equilibrium temperature is a weighted average
            T = sum(obj.quantity .* speciesTemperatures .* CpOrCv_0) / sum(obj.quantity .* CpOrCv_0);
            
            % Initialize iteration counter and convergence check variable
            it = 0;
            STOP = 1.0;
            
            % Newton-Raphson iterative loop
            while STOP > tol0 && it < itMax
                % Update iteration
                it = it + 1;
                
                % Evaluate the enthalpy/internal-energy at the current guess temperature
                HorE = obj.getPropertyListSpecies(funHorE, T);
                
                % Evaluate the function f, i.e., the difference between the initial and current enthalpy/internal-energy
                f = sum(obj.quantity .* HorE_0) - sum(obj.quantity .* HorE);
                f_rel = f / sum(obj.quantity .* HorE);
                
                % Compute the derivative using the appropriate specific heat function
                CpOrCv = obj.getPropertyListSpecies(funCpOrCv, T);
                df = -sum(obj.quantity .* CpOrCv);
                
                % Compute the temperature correction
                DeltaT = -f / df;
                
                % Update temperature
                T = T + DeltaT;
                
                % Compute STOP criterion
                STOP = max(abs([DeltaT, f_rel]));
            end

        end

        function value = getPropertyListSpecies(obj, fun, temperatures)
            % This helper method evaluates a given property function (e.g., specific heat, enthalpy)
            % for each species in the mixture at the specified temperature(s).
            %
            % Args:
            %       fun (function handle): Function to compute the property (e.g., @getHeatCapacityPressure,
            %                              @getHeatCapacityVolume, @getEnthalpy, or @getInternalEnergy).
            %       temperatures (float): Temperature(s) [K] at which to evaluate the property. If provided
            %                             as a vector, each species is evaluated at its corresponding temperature.
            %
            % Returns:
            %       value (vector): Vector containing the computed property for each species.
            %
            % Example:
            %       cp_values = obj.getPropertyListSpecies(@getHeatCapacityPressure, [300, 400, 350]);
            
            % Definitions
            numSpecies = obj.numSpecies;
            
            % If a single temperature is provided, replicate it for all species.
            if isscalar(temperatures)
                temperatures = repmat(temperatures, 1, numSpecies);
            end

            % Evaluate property for each species
            for i = 1:numSpecies
                species = obj.chemicalSystem.database.species.(obj.listSpecies{i});
                value(i) = fun(species, temperatures(i));
            end

        end

    end
    
    methods (Access = private, Static)

        function mass = getMass(system, propertiesMatrix)
            % Compute mass mixture [kg]
        
            % Compute total number of moles [mol]
            N = sum(propertiesMatrix(:, system.ind_ni));
            % Compute mean molecular weight [kg/mol]
            W = dot(propertiesMatrix(:, system.ind_ni), propertiesMatrix(:, system.ind_W)) / N;
            % Compute mass of the mixture [kg]
            mass = W * N;
        end
    
    end

    methods (Hidden)

        function obj = setPropertiesMatrixFast(obj, listSpecies, quantity, index, h0)
            % Set species and quantity and compute thermodynamic properties
           
            % Update local listSpecies and local quantity
            % obj.listSpecies = listSpecies;
            % obj.quantity = quantity;

            % Assign values to the propertiesMatrix
            setPropertiesMatrix(obj.chemicalSystem, listSpecies, quantity, obj.T, index, h0);

            % Compute thermodynamic properties
            computeProperties(obj);
        end

        function FLAG_FIXED = checkTemperatureSpecies(obj)
            % Check if condensed species can be only evaluated a particular temperature

            % Initialization
            FLAG_FIXED = false;

            for i = 1:obj.numSpecies
                % Get Species object
                species = obj.chemicalSystem.database.species.(obj.listSpecies{i});
                
                % Check if condensed species can be only evaluated a particular temperature
                if isscalar(species.T)
                    obj.Tspecies(i) = obj.chemicalSystem.database.species.(obj.listSpecies{i}).T;
                    FLAG_FIXED = true;
                end
        
            end
        
        end

    end

end