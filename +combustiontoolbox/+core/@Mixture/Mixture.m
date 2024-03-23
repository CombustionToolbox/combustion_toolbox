classdef Mixture

    % properties
    %     T (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     p (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     N (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     hf (1,1) double {mustBeNumeric}
    %     ef (1,1) double {mustBeNumeric}
    %     h (1,1) double {mustBeNumeric}
    %     e (1,1) double {mustBeNumeric}
    %     g (1,1) double {mustBeNumeric}
    %     s (1,1) double {mustBeNumeric}
    %     cp (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     cv (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     gamma (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     gamma_s (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     sound (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     s0 (1,1) double {mustBeNumeric}
    %     DhT (1,1) double {mustBeNumeric}
    %     DeT (1,1) double {mustBeNumeric}
    %     Ds (1,1) double {mustBeNumeric}
    %     rho (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     v (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     W (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     MW (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     mi (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     Xi (:, 1) double {mustBeNumeric, mustBeNonnegative}
    %     Yi (:, 1) double {mustBeNumeric, mustBeNonnegative}
    %     phase (:, 1) logical 
    %     dVdT_p (1,1) double {mustBeNumeric}
    %     dVdp_T (1,1) double {mustBeNumeric}
    %     equivalenceRatio (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     natomElements (1,:) double {mustBeNumeric, mustBeNonnegative}
    %     natomElementsReact (1,:) double {mustBeNumeric, mustBeNonnegative}
    %     errorMoles (1, 1) double {mustBeNumeric}
    %     errorMolesIons (1, 1) double {mustBeNumeric}
    %     errorProblem (1, 1) double {mustBeNumeric}
    %     chemicalSystem combustiontoolbox.core.ChemicalSystem
    %     rangeName (1, :) char {mustBeText}
    %     rangeValue (1, :) double {mustBeNumeric}
    % end

    properties
        T (1,1)
        p (1,1)
        N
        hf
        ef
        h
        e
        g
        s
        cp
        cv
        gamma
        gamma_s
        sound
        s0
        DhT
        DeT
        Ds
        rho
        v
        W
        MW
        mi
        Xi 
        Yi 
        phase
        dVdT_p
        dVdp_T
        equivalenceRatio
        natomElements
        natomElementsReact
        errorMoles
        errorMolesIons
        errorProblem
        chemicalSystem
        rangeName
        rangeValue
    end

    properties (Access = private)
        quantity
        listSpecies
        listSpeciesFuel
        listSpeciesOxidizer
        listSpeciesInert
        molesFuel
        molesOxidizer
        molesInert
        FLAGS_PROP
        pressureEOS = @eos_ideal_p; % Equation of State to compute pressure [Pa]
        volumeEOS = @eos_ideal;     % Equation of State to compute molar volume [m3/mol]
        chemicalPotentialImpEOS = @mu_imp_ideal; % Compute non ideal contribution of the chemical potential (depends of the Equation of State) [J/mol]
    end

    % properties (Hidden)
    %     cp_f (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     cp_r (1,1) double {mustBeNumeric, mustBeNonnegative}
    %     dNi_T (:, 1) double {mustBeNumeric}
    %     dN_T (1,1) double {mustBeNumeric}
    %     dNi_p (:, 1) double {mustBeNumeric}
    %     dN_p (1,1) double {mustBeNumeric}
    %     FLAG_REACTION = false
    % end
    
    properties (Hidden)
        cp_f
        cp_r
        dNi_T
        dN_T
        dNi_p
        dN_p
        FLAG_REACTION = false
    end

    methods

        function obj = Mixture(chemicalSystem, varargin)
            % The chemical system has the database, we only need the properties matrix to perform all the calculations
            % Now to be compatible with all the capabilities of the previous version, we will have to be able to add a
            % list of species for calculations and another list with all the species involved (reactants and products)

            % Definitions
            defaultTemperature = 300; % [K]
            defaultPressure = 1;      % [bar]

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'chemicalSystem'); % @(x) isa(x, 'combustiontoolbox.core.ChemicalSystem ')
            addOptional(ip, 'T', defaultTemperature, @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'p', defaultPressure, @(x) isnumeric(x) && x >= 0);
            parse(ip, chemicalSystem, varargin{:});

            % Assign properties matrix
            obj.chemicalSystem = ip.Results.chemicalSystem;

            % Get indexListSpeciesOriginal
            obj.chemicalSystem.indexListSpeciesOriginal = find_ind(obj.chemicalSystem.listSpecies, obj.chemicalSystem.listSpeciesOriginal);
            
            % Get indexReact
            obj.chemicalSystem = obj.chemicalSystem.set_react_index(obj.listSpeciesInert);

            % Compute thermodynamic properties
            % obj = obj.compute_properties(obj.chemicalSystem, ip.Results.T, ip.Results.p);
        end

        function obj = setTemperature(obj, T, varargin)
            % Compute thermodynamic properties at temperature T [K]
            
            % Definitions
            defaultUnits = 'K';

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'T', @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'units', defaultUnits, @(x) ischar(x));
            parse(ip, T, varargin{:});
            
            % Assign temperature
            obj.T = T;

            % Assign values to the propertiesMatrix
            obj.chemicalSystem = obj.chemicalSystem.set_propertiesMatrix(obj.listSpecies, obj.quantity, obj.T);
            
            % Check if initial state is defined (temperature, pressure, and composition)
            if ~sum(obj.quantity) || ~obj.p
                return
            end
            
            % Compute thermodynamic properties
            obj = obj.compute_properties(obj.chemicalSystem, T, obj.p);
        end

        function obj = setPressure(obj, p, varargin)
            % Compute thermodynamic properties at pressure p [K]
            
            % Definitions
            defaultUnits = 'bar';

            % Parse inputs
            ip = inputParser;
            addRequired(ip, 'T', @(x) isnumeric(x) && x >= 0);
            addOptional(ip, 'units', defaultUnits, @(x) ischar(x));
            parse(ip, p, varargin{:});
            
            % Assign pressure
            obj.p = p;

            % Assign values to the propertiesMatrix
            obj.chemicalSystem = obj.chemicalSystem.set_propertiesMatrix(obj.listSpecies, obj.quantity, obj.T);

            % Check if initial state is defined (temperature, pressure, and composition)
            if ~sum(obj.quantity) || ~obj.T
                return
            end

            % Compute thermodynamic properties
            obj = obj.compute_properties(obj.chemicalSystem, obj.T, p);
        end

        function obj = set(obj, listSpecies, varargin)
            % 
            
            % Definitions
            quantityDefault = 1;
            unitsDefault = 'mol';

            % 
            if nargin > 3
                type = varargin{1};
                quantity = varargin{2};

                switch lower(type)
                    case 'fuel'
                        obj.listSpeciesFuel = [obj.listSpeciesFuel, listSpecies];
                        obj.molesFuel = [obj.molesFuel, quantity];
                    case 'oxidizer'
                        obj.listSpeciesOxidizer = [obj.listSpeciesOxidizer, listSpecies];
                        obj.molesOxidizer = [obj.molesOxidizer, quantity];
                    case 'inert'
                        obj.listSpeciesInert = [obj.listSpeciesInert, listSpecies];
                        obj.molesInert = [obj.molesInert, quantity];
                end

            elseif nargin > 2
                quantity = varargin{1};
            end

            % Update local listSpecies and local quantity
            obj.listSpecies = [obj.listSpecies, listSpecies];
            obj.quantity = [obj.quantity, quantity];

            % Assign values to the propertiesMatrix
            obj.chemicalSystem = obj.chemicalSystem.set_propertiesMatrix(listSpecies, quantity, obj.T);
            
            % Check if initial state is defined (temperature, pressure, and composition)
            if ~obj.T || ~obj.p
                return
            end

            % Compute thermodynamic properties
            obj = obj.compute_properties(obj.chemicalSystem, obj.T, obj.p);

        end
        
        function obj = set_prop(obj, varargin)
            % Assign property values to the respective variables
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %
            % Optional Args:
            %     * field (str): Fieldname in Problem Description (PD)
            %     * value (float): Value/s to assing in the field in Problem Description (PD)
            %
            % Returns:
            %     self (struct): Data of the mixture, conditions, and databases
        
            try
                for i = 1:2:nargin - 1
                    self = set_prop_common(self, varargin{i}, varargin{i + 1});
                end
        
            catch
                error('Error number inputs @set_prop');
            end
        
        end
        
        % SUB-PASS FUNCTIONS
        function obj = set_prop_common(obj, name, value)
            % Assign property values to the given variable name
        
            % If the value is a char, convert it to a float - (GUI)
            if value(1) == '['
                value = sscanf(value, '[%f:%f:%f]');
                value = value(1):value(2):value(3);
            elseif sum(value == ':') == 2
                value = sscanf(value, '%f:%f:%f');
                value = value(1):value(2):value(3);
            elseif ischar(value)
                value = sscanf(value, '%f');
            end
        
            % Define flags
            if length(value) > 1
                obj.FLAGS_PROP.(name) = true;
            else
                obj.FLAGS_PROP.(name) = false;
            end
        
            % Set value
            obj.(name).value = value;
        end

        function print(obj)
            % 
            print_mixture(obj);
        end
        
    end

    methods (Access = private)

        function obj = compute_properties(obj, system, T, p)
            % Compute properties from the given chemicalSystem at pressure p [bar]
            % and temperature T [K]
            %
            % Args:
            %     obj (Mixture): Mixture class
            %     system (ChemicalSystem): 
            %     temperature (float): Temperature [K]
            %     pressure (float): Pressure [bar]
            %
            % Returns:
            %     obj (Mixture): Mixture class with the computed properties
            %
            % Example:
            %     mix = compute_properties(mix, system, T, p)
        
            % Definitions
            R0 = 8.31446261815324; % [J/(K mol)]. Universal gas constant
            propertiesMatrix = system.propertiesMatrix; % Properties matrix

            % Initialization
            obj.errorMoles = 0;
            obj.errorMolesIons = 0;

            % Inputs
            obj.T = T; % [K]
            obj.p = p; % [bar]

            % Unpack propertiesMatrix
            Ni = propertiesMatrix(:, system.ind_ni); % [mol]
            obj.N = sum(propertiesMatrix(:, system.ind_ni)); % [mol]
            obj.hf = dot(propertiesMatrix(:, system.ind_hfi), Ni); % [kJ]
            obj.h = dot(propertiesMatrix(:, system.ind_hi), Ni); % [kJ]
            obj.ef = dot(propertiesMatrix(:, system.ind_efi), Ni); % [kJ]
            obj.cp = dot(propertiesMatrix(:, system.ind_cpi), Ni); % [J/K]
            obj.s0 = dot(propertiesMatrix(:, system.ind_si), Ni); % [kJ/K]
            obj.phase = propertiesMatrix(:, system.ind_phase); % [bool]

            % Compute total composition of gas species [mol]
            N_gas = sum(Ni(obj.phase == 0));

            % Compute molar fractions [-]
            obj.Xi = Ni / obj.N;

            % Compute molecular weight (gases) [g/mol]
            obj.W = dot(Ni, propertiesMatrix(:, system.ind_W)) / N_gas;

            % Compute mean molecular weight [g/mol]
            obj.MW = dot(Ni, propertiesMatrix(:, system.ind_W)) / obj.N;

            % Compute mass mixture [kg]
            obj.mi = obj.MW * obj.N * 1e-3; % [kg]

            % Compute mass fractions [-]
            obj.Yi = (Ni .* propertiesMatrix(:, system.ind_W) * 1e-3) ./ obj.mi;

            % Get non zero species
            FLAG_NONZERO = obj.Xi > 0;

            % Compute vector atoms of each element
            obj.natomElements = sum(Ni .* system.stoichiometricMatrix, 1);

            % Compute vector atoms of each element without frozen species
            obj.natomElementsReact = sum(propertiesMatrix(system.indexReact, system.ind_ni) .* system.stoichiometricMatrix(system.indexReact, :), 1);

            % Compute volume [m3]
            obj.v = obj.volumeEOS(T, convert_bar_to_Pa(p), obj.chemicalSystem.listSpecies, obj.Xi) * N_gas;

            % Compute density [kg/m3]
            obj.rho = obj.mi / obj.v;

            % Compute internal energy [kJ]
            obj.e = obj.h - N_gas * R0 * T * 1e-3;

            % Compute thermal internal energy [kJ]
            obj.DeT = obj.e - obj.ef;

            % Compute thermal enthalpy [kJ]
            obj.DhT = obj.h - obj.hf;

            % Compute entropy of mixing [kJ/K]
            obj.Ds = obj.compute_entropy_mixing(Ni, N_gas, R0, FLAG_NONZERO);

            % Compute entropy [kJ/K]
            obj.s = obj.s0 + obj.Ds;

            % Compute Gibbs energy [kJ]
            obj.g = obj.h - obj.T * obj.s;

            % Compute specific heat at constant volume [J/K]
            obj.cv = obj.cp - R0 * N_gas;

            % Compute Adibatic index [-]
            obj.gamma = obj.cp / obj.cv;
            
            % Compute sound velocity [m/s]
            obj.sound = sqrt(obj.gamma * convert_bar_to_Pa(p) / obj.rho);
            
            % Correction of: cP, cv, gamma, and speed of sound as consequence of the
            % chemical reaction
            if obj.FLAG_REACTION % isfield(self, 'dNi_T')
                obj.dVdT_p = obj.dN_T + 1; % [-]
                obj.dVdp_T = obj.dN_p - 1; % [-]
        
                if ~any(isnan(obj.dNi_T)) && ~any(isinf(obj.dNi_T))
                    delta = ~obj.phase;
                    h0_j = propertiesMatrix(:, system.ind_hi) * 1e3; % [J/mol]
                    obj.cp_r = sum(h0_j / T .* (1 + delta .* (Ni - 1)) .* obj.dNi_T, 'omitnan'); % [J/K]
                    obj.cp_f = obj.cp;
                    obj.cp = obj.cp_f + obj.cp_r; % [J/K]
                    obj.cv = obj.cp + (N_gas * R0 * obj.dVdT_p^2) / obj.dVdp_T; % [J/K]
                    obj.gamma = obj.cp / obj.cv; % [-]
                    obj.gamma_s =- obj.gamma / obj.dVdp_T; % [-]
                    obj.sound = sqrt(obj.gamma_s * convert_bar_to_Pa(p) / obj.rho); % [m/s]

                    return
                end

                obj.gamma_s =- 1 / obj.dVdp_T; % [-]
                
                return
            end
            
            obj.dVdT_p = 1;          % [-]
            obj.dVdp_T = -1;         % [-]
            obj.gamma_s = obj.gamma; % [-]
        
        end

        function Ds = compute_entropy_mixing(obj, Ni, N_gas, R0, FLAG_NONZERO)
            % Compute entropy of mixing [kJ/K]
            %
            % Args:
            %     obj (Mixture): Mixture class
            %     Ni (moles): Vector of moles of each species [mol]
            %     N_gas (moles): Total moles of gas species [mol]
            %     R0 (float): Universal gas constant [J/(mol-K)]
            %     FLAG_NONZERO (bool): Vector of nonzero species
            %
            % Returns:
            %     DS (float): Entropy of mixing [kJ/K]
            %
            % Note: 
            %     only nonzero for gaseous species
            %
            % Example:
            %     Ds = compute_entropy_mixing(mix, Ni, N_gas, R0, FLAG_NONZERO)
        
            Dsi = Ni(FLAG_NONZERO) .* log(Ni(FLAG_NONZERO) / N_gas * obj.p) .* (1 - obj.phase(FLAG_NONZERO));
            Ds = -R0 * sum(Dsi) * 1e-3;
        end

    end

    methods (Hidden)

        function obj = set_fast(obj, listSpecies, quantity, index, h0)
            % 
           
            % Update local listSpecies and local quantity
            obj.listSpecies = listSpecies;
            obj.quantity = quantity;

            % Assign values to the propertiesMatrix
            obj.chemicalSystem = obj.chemicalSystem.set_propertiesMatrix(listSpecies, quantity, obj.T, index, h0);

            % Compute thermodynamic properties
            obj = obj.compute_properties(obj.chemicalSystem, obj.T, obj.p);

        end

    end

end