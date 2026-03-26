classdef EquationStatePengRobinson < combustiontoolbox.core.EquationState
    % The :mat:class:`EquationStatePengRobinson` class implements the 
    % Peng-Robinson equation of state for real gases.
    %
    % Example:
    %      eos = EquationStatePengRobinson();
    %
    % See also: :mat:class:`EquationState`, :mat:class:`Mixture`

    properties (Access = public)
        tol0 = 1e-8;        % Tolerance for root finding
    end

    properties (Access = private)
        cachedListSpecies   % Cell array of strings to validate the cache
        temperatureCritical % Critical temperatures of all species [K]
        pressureCritical    % Critical pressures of all species [Pa]
        acentricFactor      % Acentric factors of all species [-]
        FLAG_VALID          % Flag array identifying species with valid PR data
    end

    properties (Constant, Access = private)
        R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
    end

    methods (Access = public)

        function pressure = getPressure(obj, temperature, molarVolume, molarFractions, chemicalSystem, varargin)
            % Compute pressure [Pa] using the Peng-Robinson equation of state, namely:
            %
            % .. math::
            % 
            %       `P = \\frac{RT}{V - b} - \\frac{a}{V^2 + 2bV - b^2},`
            %
            % where :math:`a` and :math:`b` are mixture parameters computed using 
            % van der Waals one-fluid mixing rules.
            %
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     pressure (float): Pressure of the mixture [Pa]
            %
            % Example:
            %     P = getPressure(obj, 300, 0.024, [0.5, 0.5], chemicalSystem)

            % Compute mixture parameters
            [a_mix, b_mix, ~, ~] = obj.getMixtureParameters(temperature, molarFractions, chemicalSystem);
            
            % Compute pressure [Pa]
            pressure = (obj.R0 * temperature) / (molarVolume - b_mix) - ...
                       a_mix / (molarVolume^2 + 2 * b_mix * molarVolume - b_mix^2);
        end
        
        function molarVolume = getVolume(obj, temperature, pressure, molarFractions, chemicalSystem, varargin)
            % Compute gas-phase molar volume [m3/mol] by solving the Peng-Robinson
            % cubic equation of state for the given temperature and pressure
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %
            % Example:
            %     V = getVolume(obj, 300, 1e5, [0.5, 0.5], chemicalSystem)

            % Compute mixture parameters
            [a_mix, b_mix, ~, ~] = obj.getMixtureParameters(temperature, molarFractions, chemicalSystem);
            
            % Dimensionless coefficients A and B
            A = (a_mix * pressure) / (obj.R0^2 * temperature^2);
            B = (b_mix * pressure) / (obj.R0 * temperature);
            
            % Cubic coefficients for Z^3 + c2*Z^2 + c1*Z + c0 = 0
            coeffs = [1.0, -(1.0 - B), (A - 2*B - 3*B^2), -(A*B - B^2 - B^3)];
            
            % Solve for Z and pick the largest real root (gas phase)
            Z_roots = roots(coeffs);
            Z_real = real(Z_roots(abs(imag(Z_roots)) < obj.tol0));
            
            if isempty(Z_real)
                error('EquationStatePengRobinson:getVolume', 'No real roots found for Z.');
            end
            
            Z_gas = max(Z_real);

            % Compute molar volume [m3/mol]
            molarVolume = (Z_gas * obj.R0 * temperature) / pressure;
        end

        function [dPdV_T, dPdT_V] = getPressureDerivativesDimensional(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
            % Compute dimensional partial pressure derivatives for the mixture assuming frozen chemistry
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * dPdV_T (float): Partial derivative of pressure with respect to volume at constant temperature [Pa/(m3/mol)]
            %     * dPdT_V (float): Partial derivative of pressure with respect to temperature at constant volume [Pa/K]
            %
            % Example:
            %     [dPdV_T, dPdT_V] = getPressureDerivativesDimensional(obj, 300, 1e5, 0.024, [0.5, 0.5], chemicalSystem)

            % Compute mixture parameters
            [a_mix, b_mix, dadT_mix, ~] = obj.getMixtureParameters(temperature, molarFractions, chemicalSystem);
            
            % Compute dimensional pressure derivatives
            dPdV_T = -(obj.R0 * temperature) / (molarVolume - b_mix)^2 + (2 * a_mix * (molarVolume + b_mix)) / (molarVolume^2 + 2*b_mix*molarVolume - b_mix^2)^2;
            dPdT_V = obj.R0 / (molarVolume - b_mix) - dadT_mix / (molarVolume^2 + 2*b_mix*molarVolume - b_mix^2);
        end

        function [heatCapacityPressureDeparture, enthalpyDeparture, entropyDeparture] = getDepartureFunctions(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
            % Compute thermodynamic departure functions for the mixture using the Peng-Robinson equation of state
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * heatCapacityPressureDeparture (float): Heat capacity at constant pressure departure [J/(mol-K)]
            %     * enthalpyDeparture (float): Enthalpy departure [J/mol]
            %     * entropyDeparture (float): Entropy departure [J/(mol-K)]
            %
            % Example:
            %     [dcp, dh, ds] = getDepartureFunctions(obj, 300, 1e5, 0.024, [0.5, 0.5], chemicalSystem)

            % Compute mixture parameters
            [a_mix, b_mix, dadT_mix, d2adT2_mix] = obj.getMixtureParameters(temperature, molarFractions, chemicalSystem);
            
            % If mixture behaves ideally (b_mix is zero), return zeros for all departure functions
            if b_mix < 1e-15
                heatCapacityPressureDeparture = 0;
                enthalpyDeparture = 0; 
                entropyDeparture = 0; 
                return
            end

            Z = obj.getCompressibilityFactor(temperature, pressure, molarVolume);
            B = (b_mix * pressure) / (obj.R0 * temperature);
            
            % Common terms for departure functions
            arg = (Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B);
            logTerm = log(max(arg, 1e-12));
            denom = 2 * sqrt(2) * b_mix;

            % Enthalpy departure [J/mol]
            enthalpyDeparture = obj.R0 * temperature * (Z - 1) + ((temperature * dadT_mix - a_mix) / denom) * logTerm;
            
            % Entropy departure [J/(mol-K)]
            entropyDeparture = obj.R0 * log(max(Z - B, 1e-12)) + (dadT_mix / denom) * logTerm;
            
            % Heat capacity at constant volume departure [J/(mol-K)]
            heatCapacityVolumeDeparture = (temperature * d2adT2_mix / denom) * logTerm;

            % Compute pressure derivatives
            [dPdV_T, dPdT_V] = obj.getPressureDerivativesDimensional(temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin{:});
            
            % Heat capacity at constant pressure departure [J/(mol-K)]
            heatCapacityPressureDeparture = heatCapacityVolumeDeparture + (-temperature * (dPdT_V^2) / dPdV_T) - obj.R0;
        end

    end
    
    methods (Access = private)
        
        function initializeCache(obj, chemicalSystem)
            % Cache the database values once to avoid field lookups during iterative solver loops
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Example:
            %     obj.initializeCache(chemicalSystem)
            
            % Definitions
            listSpecies = chemicalSystem.listSpecies;
            numSpecies = chemicalSystem.numSpecies;
            
            % Preallocate arrays
            Tc = zeros(1, numSpecies);
            Pc = zeros(1, numSpecies);
            omega = zeros(1, numSpecies);
            
            % Extract Species objects from the chemical system
            species = chemicalSystem.species;
            for i = 1:numSpecies
                name = listSpecies{i};
                Tc(i) = species.(name).Tcritical;
                Pc(i) = species.(name).Pcritical;
                omega(i) = species.(name).acentricFactor;
            end
            
            obj.temperatureCritical = Tc;    % [K]
            obj.pressureCritical = Pc * 1e5; % [Pa]
            obj.acentricFactor = omega;      % [-]
            
            % Identify species with valid PR data (not NaN and > 0)
            obj.FLAG_VALID = ~isnan(Tc) & (Tc > 0) & ~isnan(Pc) & (Pc > 0);
            
            % Cache the list of species to validate future calls
            obj.cachedListSpecies = listSpecies; 
        end

        function [a_mix, b_mix, dadT_mix, d2adT2_mix] = getMixtureParameters(obj, temperature, molarFractions, chemicalSystem)
            % Compute mixture parameters using van der Waals one-fluid mixing rules
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * a_mix (float): Mixture attraction parameter [J-m3/mol^2]
            %     * b_mix (float): Mixture co-volume parameter [m3/mol]
            %     * dadT_mix (float): First temperature derivative of a_mix [J-m3/(mol^2 K)]
            %     * d2adT2_mix (float): Second temperature derivative of a_mix [J-m3/(mol^2 K^2)]
            %
            % Example:
            %     [a_mix, b_mix, dadT_mix, d2adT2_mix] = getMixtureParameters(obj, 300, [0.5, 0.5], chemicalSystem)
            
            % Rebuild cache if empty or if species list changed/reordered
            if isempty(obj.cachedListSpecies) || ~isequal(obj.cachedListSpecies, chemicalSystem.listSpecies)
                obj.initializeCache(chemicalSystem);
            end
            
            % Find species that are both active (>0) AND have valid PR data
            FLAG_ACTIVE = (molarFractions(:) > 0) & obj.FLAG_VALID(:);
            
            % If no real species are present, mixture is purely ideal
            if ~any(FLAG_ACTIVE)
                a_mix = 0; b_mix = 0; dadT_mix = 0; d2adT2_mix = 0;
                return;
            end
            
            % Extract data for active species
            X_active = molarFractions(FLAG_ACTIVE);
            Tc = obj.temperatureCritical(FLAG_ACTIVE);
            Pc = obj.pressureCritical(FLAG_ACTIVE);
            omega = obj.acentricFactor(FLAG_ACTIVE);
            
            % Pure species parameters
            kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega.^2;
            Tr = temperature ./ Tc;
            sqrtTr = sqrt(Tr);
            alpha = (1 + kappa .* (1 - sqrtTr)).^2;
            
            a_i = 0.45724 * (obj.R0^2 * Tc.^2 ./ Pc) .* alpha;
            b_i = 0.07780 * (obj.R0 * Tc ./ Pc);
            
            % Temperature derivatives
            dalpha_dT = -kappa ./ (sqrtTr .* Tc) .* (1 + kappa .* (1 - sqrtTr));
            d2alpha_dT2 = kappa .* (kappa + 1) ./ (2 * Tc.^2 .* Tr.^(3/2));
            
            a0 = 0.45724 * (obj.R0^2 * Tc.^2 ./ Pc);
            da_dT_i = a0 .* dalpha_dT;
            d2a_dT2_i = a0 .* d2alpha_dT2;
            
            % Apply mixing rules
            b_mix = dot(X_active, b_i);
            
            sqrt_a_i = sqrt(max(a_i, 1e-20));
            S1 = dot(X_active, sqrt_a_i);
            a_mix = S1^2;
            
            S2 = dot(X_active, da_dT_i ./ (2 * sqrt_a_i));
            dadT_mix = 2 * S1 * S2;
            
            S3 = dot(X_active, (d2a_dT2_i ./ (2 * sqrt_a_i)) - (da_dT_i.^2 ./ (4 * sqrt_a_i.^3)));
            d2adT2_mix = 2 * S2^2 + 2 * S1 * S3;
        end

        function [temperatureCritical_mix, pressureCritical_mix, acentricFactor_mix] = getPseudoCriticalProperties(obj, molarFractions, chemicalSystem)
            % Computes pseudo-critical properties for multi-component mixtures
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * temperatureCritical_mix (float): Pseudo-critical temperature of the mixture [K]
            %     * pressureCritical_mix (float): Pseudo-critical pressure of the mixture [Pa]
            %     * acentricFactor_mix (float): Pseudo-critical acentric factor of the mixture [-]
            %
            % Example:
            %     [temperatureCritical_mix, pressureCritical_mix, acentricFactor_mix] = getPseudoCriticalProperties(obj, [0.5, 0.5], chemicalSystem)
            
            if isempty(obj.cachedListSpecies) || ~isequal(obj.cachedListSpecies, chemicalSystem.listSpecies)
                obj.initializeCache(chemicalSystem);
            end
            
            % Definitions
            mask = (molarFractions(:) > 0) & obj.FLAG_VALID(:);
            Xi = molarFractions(mask);
            Tc_i = obj.temperatureCritical(mask);
            Pc_i = obj.pressureCritical(mask);
            omega_i = obj.acentricFactor(mask);
            
            % a and b at critical point (alpha = 1)
            a_i_tc = 0.45724 * (obj.R0^2 * Tc_i.^2 ./ Pc_i);
            b_i = 0.07780 * (obj.R0 * Tc_i ./ Pc_i);
            
            % Mixture parameters at critical condition
            a_mix_tc = ( dot(Xi, sqrt(a_i_tc)) )^2;
            b_mix = dot(Xi, b_i);
            
            % Back-calculate pseudo-critical T and P
            temperatureCritical_mix = (a_mix_tc * 0.07780) / (b_mix * 0.45724 * obj.R0);
            pressureCritical_mix = (0.07780 * obj.R0 * temperatureCritical_mix) / b_mix;
            acentricFactor_mix = dot(Xi, omega_i);
        end
        
    end
end