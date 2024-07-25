classdef EquilibriumSolver < handle
    % EquilibriumSolver: A class for solving chemical equilibrium problems
    
    properties
        problemType                % Problem type [TP, TV, HP, EV, SP, SV]
        % * Chemical equilibrium TP, TV (CT-EQUIL module)
        tolGibbs = 1e-6            % Tolerance of the Gibbs/Helmholtz minimization method
        tolE = 1e-6                % Tolerance of the mass balance
        tolMoles = 1e-14           % Tolerance of the composition of the mixture                       
        tolMolesGuess = 1e-6       % Tolerance of the molar composition of the mixture (guess)
        tolMultiplierIons = 1e-4   % Tolerance of the dimensionless Lagrangian multiplier - ions
        tolTau = 1e-25             % Tolerance of the slack variables for condensed species
        itMaxGibbs = 70            % Max number of iterations - Gibbs/Helmholtz minimization method
        itMaxIons = 30             % Max number of iterations - charge balance (ions)
        temperatureIons = 0        % Minimum temperature [K] to consider ionized species
        % * Chemical equilibrium HP, EV, SP, SV (CT-EQUIL module)
        tol0 = 1e-3                % Tolerance of the root finding algorithm
        itMax = 30                 % Max number of iterations - root finding method
        rootMethod = @newton       % Root finding method [newton (2nd order), steff (2nd order), or nsteff (3rd order)]
        root_T0_l = 1000           % First temperature guess [K] left branch - root finding method
        root_T0_r = 3000           % First temperature guess [K] right branch - root finding method
        root_T0   = 3000           % Temperature guess [K] if it's outside previous range - root finding method
        % * Flags
        FLAG_EXTRAPOLATE = true    % Flag indicating linear extrapolation of the polynomials fits
        FLAG_FAST = true           % Flag indicating use guess composition of the previous computation
        FLAG_TCHEM_FROZEN = false  % Flag to consider a thermochemically frozen gas (calorically perfect gas)
        FLAG_FROZEN = false        % Flag to consider a calorically imperfect gas with frozen chemistry
        FLAG_EOS = false           % Flag to use non-ideal Equation of States (EoS)
        FLAG_RESULTS = true        % Flag to print results
        FLAG_TIME = true           % Flag to print elapsed time
        % * Miscellaneous
        time
    end

    methods

        function obj = EquilibriumSolver(varargin)
            % Constructor
            defaultProblemType = 'TP';

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'TP', 'TV', 'HP', 'EV', 'SP', 'SV'})));
            addParameter(p, 'tolGibbs', obj.tolGibbs, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'tolE', obj.tolE, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'tolMoles', obj.tolMoles, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'tolMolesGuess', obj.tolMolesGuess, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'tolMultiplierIons', obj.tolMultiplierIons, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'tolTau', obj.tolTau, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'itMaxGibbs', obj.itMaxGibbs, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'itMaxIons', obj.itMaxIons, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'temperatureIons', obj.temperatureIons, @(x) isnumeric(x) && x >= 0);
            addParameter(p, 'tol0', obj.tol0, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'itMax', obj.itMax, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'rootMethod', obj.rootMethod, @(method) ismember(method, {@newton, @steff, @nsteff}));
            addParameter(p, 'root_T0_l', obj.root_T0_l, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'root_T0_r', obj.root_T0_r, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'root_T0', obj.root_T0, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'FLAG_EXTRAPOLATE', obj.FLAG_EXTRAPOLATE, @(x) islogical(x));
            addParameter(p, 'FLAG_FAST', obj.FLAG_FAST, @(x) islogical(x));
            addParameter(p, 'FLAG_TCHEM_FROZEN', obj.FLAG_TCHEM_FROZEN, @(x) islogical(x));
            addParameter(p, 'FLAG_FROZEN', obj.FLAG_FROZEN, @(x) islogical(x));
            addParameter(p, 'FLAG_EOS', obj.FLAG_EOS, @(x) islogical(x));
            addParameter(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.tolGibbs = p.Results.tolGibbs;
            obj.tolE = p.Results.tolE;
            obj.tolMoles = p.Results.tolMoles;
            obj.tolMolesGuess = p.Results.tolMolesGuess;
            obj.tolMultiplierIons = p.Results.tolMultiplierIons;
            obj.tolTau = p.Results.tolTau;
            obj.itMaxGibbs = p.Results.itMaxGibbs;
            obj.itMaxIons = p.Results.itMaxIons;
            obj.temperatureIons = p.Results.temperatureIons;
            obj.tol0 = p.Results.tol0;
            obj.itMax = p.Results.itMax;
            obj.rootMethod = p.Results.rootMethod;
            obj.root_T0_l = p.Results.root_T0_l;
            obj.root_T0_r = p.Results.root_T0_r;
            obj.root_T0 = p.Results.root_T0;
            obj.FLAG_EXTRAPOLATE = p.Results.FLAG_EXTRAPOLATE;
            obj.FLAG_FAST = p.Results.FLAG_FAST;
            obj.FLAG_TCHEM_FROZEN = p.Results.FLAG_TCHEM_FROZEN;
            obj.FLAG_FROZEN = p.Results.FLAG_FROZEN;
            obj.FLAG_EOS = p.Results.FLAG_EOS;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
        end

        function mix = solve(obj, mix, varargin)
            % Obtain chemical equilibrium composition and thermodynamic properties
            if nargin > 2
                obj.equilibrate(mix, varargin{1});
            else
                obj.equilibrate(mix);
            end
            
            % Set problemType
            mix.problemType = obj.problemType;
            
            % Print results
            if obj.FLAG_RESULTS
                print(mix);
            end

        end

        function mixArray = solveArray(obj, mixArray, varargin)
            % Obtain chemical equilibrium composition and thermodynamic properties for an array of values
            
            % Definitions
            n = length(mixArray);
            
            % Timer
            obj.time = tic;

            % Calculations
            obj.solve(mixArray(n));
            
            for i = n-1:-1:1
                obj.solve(mixArray(i), mixArray(i + 1));
            end

            % Timer
            obj.time = toc(obj.time);

            % Print elapsed time
            printTime(obj);
        end

        function mix2 = equilibrate(obj, mix2, varargin)
            % Obtain properties at equilibrium for the given thermochemical transformation
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix (Mixture): Mixture considering a thermochemical frozen gas
            %
            % Returns:
            %     mix (Mixture): Mixture at chemical equilibrium for the given thermochemical transformation
            %
            % Example:
            %     * mix = equilibrate(EquilibriumSolver(), 'TP', mix)
            
            % Initialization
            mix1 = mix2.copy();

            if obj.FLAG_TCHEM_FROZEN
                % TO BE IMPLEMENTED
                mix2.errorProblem = 0;
                return
            end
            
            % Default
            mix_guess = [];

            % Unpack guess
            if nargin > 2
                mix_guess = varargin{1};
            end

            % Get attribute xx of the specified transformations
            attr_name = get_attr_name(obj);
            % Compute initial guess
            [guess, guess_moles] = get_guess(obj, mix1, mix2, mix_guess, attr_name);
            % If the problem type is SV, the product's volume is based on the given v_P/v_R ratio
            % mix2 = set_volume_SV(obj, mix2);
            % Root finding: find the value x that satisfies f(x) = mix2.xx(x) - mix1.xx = 0
            [T, STOP, guess_moles] = rootFinding(obj, mix1, mix2, attr_name, guess, guess_moles);
            % Compute properties
            obj.equilibrate_T(mix1, mix2, T, guess_moles);
            % Check convergence in case the problemType is TP (defined Temperature and Pressure)
            print_convergence(mix2.errorMoles, obj.tolGibbs, mix2.errorMolesIons, obj.tolMultiplierIons, obj.problemType)
            % Save error from root finding algorithm
            mix2.errorProblem = STOP;

            % SUB-PASS FUNCTIONS
            function attr_name = get_attr_name(obj)
                % Get attribute of the problem type
                switch upper(obj.problemType)
                    case {'TP', 'TV'}
                        attr_name = 'T';
                    case 'HP'
                        attr_name = 'h';
                    case 'EV'
                        attr_name = 'e';
                    case {'SP', 'SV'}
                        attr_name = 's';
                end

            end

            function [guess, guess_moles] = get_guess(obj, mix1, mix2, mix_guess, attr_name)
                % Get initial estimates for temperature and molar composition

                % Initialization
                if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                    % guess = get_transformation(obj, 'TP');
                    guess = mix2.T;
                    guess_moles = [];

                    % if ~isempty(mix_guess)
                    %     guess_moles = mix_guess.Xi * mix_guess.N;
                    % end
                    
                    return
                end

                if ~isempty(mix_guess)
                    guess = mix_guess.T;
                    guess_moles = mix_guess.Xi * mix_guess.N;
                else
                    guess = obj.regulaGuess(mix1, mix2, attr_name);
                    guess_moles = [];
                end

            end

            function [x, STOP, guess_moles] = rootFinding(obj, mix1, mix2, attr_name, x0, guess_moles)
                % Calculate the temperature value that satisfied the problem conditions
                % using the @rootMethod
                [x, STOP, guess_moles] = obj.rootMethod(obj, mix1, mix2, attr_name, x0, guess_moles);
            end

            function print_convergence(STOP, TOL, STOP_ions, TOL_ions, ProblemType)
                % Print tolerance error if the convergence criteria was not satisfied

                if ~strcmpi(ProblemType, 'TP')
                    return
                end

                if STOP > TOL
                    fprintf('***********************************************************\n')
                    fprintf('Convergence error number of moles:   %1.2e\n', STOP);
                end

                if STOP_ions > TOL_ions
                    fprintf('***********************************************************\n')
                    fprintf('Convergence error in charge balance: %1.2e\n', STOP_ions);
                end

            end

            function mix = set_volume_SV(obj, mix)
                % If the problem type is SV, the product's volume is based on the given v_P/v_R ratio
                if ~strcmpi(obj.problemType, 'SV')
                    return
                end

                %mix.v = mix.v * obj.PD.vP_vR.value;
                fprintf('\nto be clarified\n')
            end

        end

        function mix2 = equilibrate_T(obj, mix1, mix2, T, varargin)
            % Obtain equilibrium properties and composition for the given temperature [K] and pressure [bar]
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     T (float): Temperature [K]
            %
            % Optional Args:
            %     guess_moles (float): Mixture composition [mol] of a previous computation
            %
            % Returns:
            %     mix2 (Mixture): Properties of the final mixture
            %
            % Example:
            %     mix2 = equilibrate_T(EquilibriumSolver(), mix1, mix2, 3000)
            
            % Import packages
            import combustiontoolbox.utils.findIndex
            
            % Check if calculations are for a thermochemical frozen gas (calorically perfect)
            if obj.FLAG_TCHEM_FROZEN
                obj.equilibrate_T_tchem(mix1, mix2, T);
                return
            end
        
            % Check if calculations are for a calorically imperfect gas with frozen chemistry
            if obj.FLAG_FROZEN
                % Computed by default when defining the mixture (Mixture)

                % Update thermodynamic properties assuming a thermally perfect gas
                setTemperature(mix2, T);
                return
            end

            % Definitions
            N_mix0 = moles(mix1); % Get moles of inert species
            system = mix2.chemicalSystem;
            systemProducts = mix2.chemicalSystemProducts;
            
            % Unpack
            guess_moles = unpack(varargin);

            % Check flag
            if ~obj.FLAG_FAST, guess_moles = []; end

            % Compute number of moles
            [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = selectEquilibrium(obj, systemProducts, T, mix1, mix2, guess_moles);
            
            % Compute property matrix of the species at chemical equilibrium
            setMixture(mix2);

            % NESTED FUNCTIONS
            function guess_moles = unpack(value)
                % Unpack inputs
                if isempty(value)
                    guess_moles = [];
                else
                    guess_moles = value{1};
                end
            
            end
            
            function pP = computePressure(mix, T, moles, index)
                % Compute pressure [bar] of product mixture
                vMolar = vSpecific2vMolar(mix, mix.vSpecific, moles, sum(moles(system.indexGas)), index);
                pP = mix.equationOfState.getPressure(T, vMolar, system.listSpecies, mix.Xi) * 1e-5;
            end
            
            function [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = selectEquilibrium(obj, system, T, mix1, mix2, guess_moles)
                % Select equilibrium: TP: Gibbs; TV: Helmholtz
                
                if strfind(obj.problemType, 'P') == 2
                    [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumGibbs(obj, system, mix2.p, T, mix1, guess_moles);
                else
                    [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, mix2.v, T, mix1, guess_moles);
                end

            end
            
            function mix = setMixture(mix)
                % Compute properties of final mixture
                
                % Reshape composition matrix N, and partial composition partial derivatives 
                N = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, N);
                dNi_T = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, dNi_T);
                dNi_p = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, dNi_p);
                h0 = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, h0);
                % h0(system.indexFrozen) = set_h0(system.listSpecies(system.indexFrozen), T, system.species); 
                N(system.indexFrozen) = N_mix0(system.indexFrozen);

                % Assign values
                mix.T = T;
                mix.dNi_T = dNi_T; mix.dN_T = dN_T;
                mix.dNi_p = dNi_p; mix.dN_p = dN_p;
                mix.FLAG_REACTION = true;
                mix.errorMoles = STOP;
                mix.errorMolesIons = STOP_ions;

                % Clean chemical system
                system.clean;
                
                % Check if problemType is at constant volume
                if strfind(obj.problemType, 'V') == 2
                    mix.p = computePressure(mix, T, N, [system.indexGas, system.indexCondensed]);
                end
                
                % Get indexSpecies from indexProducts
                indexSpecies = findIndex(system.listSpecies, systemProducts.listSpecies(indexProducts));

                % Compute properties of final mixture
                mix.set_fast(system.listSpecies, N', [indexSpecies, system.indexFrozen], h0);
            end
            
            function vector = reshapeVector(system, index, ind_modified, vector_modified)
                % Reshape vector containing all the species
                vector = system.molesPhaseMatrix(:, 1);
                vector(index, 1) = vector_modified(ind_modified, 1);
            end
            
        end

        function mix2 = equilibrate_T_tchem(obj, mix1, mix2, T)
            % Obtain equilibrium properties and composition for the given
            % temperature [K] and pressure [bar] assuming a calorically perfect gas
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     T (float): Temperature [K]
            %
            % Returns:
            %     mix2 (Mixture): Properties of the final mixture assuming a calorically perfect gas
            %
            % Example:
            %     mix2 = equilibrate_T_tchem(EquilibriumSolver(), mix1, mix2, 3000)
        
            % Recompute properties of mix2
            setTemperature(mix2, T);

            % Change properties that remains thermochemically frozen
            mix2.cp = mix1.cp;
            mix2.cv = mix1.cv;
            mix2.gamma = mix1.gamma;
            mix2.gamma_s = mix1.gamma_s;
            mix2.sound = sqrt(mix2.gamma * convert_bar_to_Pa(mix2.p) / mix2.rho);
            
            if ~isempty(mix2.u)
                mix2.mach = mix2.u / mix2.sound;
            end
            
        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            
            if ~obj.FLAG_TIME
                return
            end

            fprintf('\nElapsed time is %.5f seconds\n', obj.time);
        end

        function plot(obj, mixArray, varargin)
            % Plot results 
            
            % Import packages
            import combustiontoolbox.utils.display.plotComposition
            import combustiontoolbox.utils.display.plotProperties
            
            % Check if is a scalar value
            if isscalar(mixArray)
                return
            end
            
            % Plot molar fractions - mixArray
            ax1 = plotComposition(mixArray(1), mixArray, mixArray(1).rangeName, 'Xi', 'mintol', 1e-14);
        
            % Plot properties - mixArray
            ax2 = plotProperties(repmat({mixArray(1).rangeName}, 1, 9), mixArray, {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []});

            for i = 1:nargin - 2
                % Unpack input
                mixArray = varargin{i};

                % Plot molar fractions - mixArray_i
                ax1 = plotComposition(mixArray(1), mixArray, mixArray(1).rangeName, 'Xi', 'mintol', 1e-14);
            
                % Plot properties - mixArray_i
                ax2 = plotProperties(repmat({mixArray(1).rangeName}, 1, 9), mixArray, {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []}, 'ax', ax2);
            end

        end

    end
    
    methods (Access = private, Static)
        [N, NP] = equilibriumGuess(N, NP, A0, muRT, b0, index, indexGas, indexIons, NG, guess_moles)
        [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N0, A0, NE, indexGas, indexCondensed, indexElements, H0RT)
        [indexCondensed, FLAG_CONDENSED, dL_dnj] = equilibriumCheckCondensed(A0, pi_i, W, indexCondensed, muRT, NC_max, FLAG_ONE, FLAG_RULE)

        function [A0, indexRemoveSpecies, ind_E, NatomE] = removeElements(NatomE, A0, ind_E, tol)
            % Find zero sum elements
        
            % Define temporal fictitious value if there are ionized species
            temp_NatomE = NatomE;
            temp_NatomE(ind_E) = 1;
        
            % Get flag of elements to be removed from stoichiometrix matrix
            FLAG_REMOVE_ELEMENTS = temp_NatomE' <= tol;
            
            % Get the species to be removed from stoichiometrix matrix
            indexRemoveSpecies = find(sum(A0(:, FLAG_REMOVE_ELEMENTS) > 0, 2) > 0);
        
            % Update stoichiometrix matrix
            A0(:, FLAG_REMOVE_ELEMENTS) = [];
        
            % Set number of atoms
            NatomE(FLAG_REMOVE_ELEMENTS) = [];
            
            % Check position "element" electron
            if ind_E
                ind_E = ind_E - sum(FLAG_REMOVE_ELEMENTS(1:ind_E-1));
            end
    
        end

        function [index, indexGas, indexCondensed, indexIons, indexElements, NE, NG, NS] = tempValues(system, NatomE)
            % List of indices with nonzero values and lengths
            indexElements = 1:length(NatomE);
            indexGas = system.indexGas;
            indexCondensed = system.indexCondensed;
            indexIons = system.indexGas(system.indexIons);
            indexCryogenic = system.indexCryogenic;
            index = [indexGas, indexCondensed];

            % Remove cryogenic species from calculations
            for i = 1:length(indexCryogenic)
                index(index == indexCryogenic(i)) = [];
                indexCondensed(indexCondensed == indexCryogenic(i)) = [];
            end
        
            % Update lengths
            NE = length(NatomE);
            NG = length(indexGas);
            NS = length(index);
        end
        
        function [index, indexCondensed, indexGas, indexIons, NG, NS, N] = updateTemp(N, index, indexCondensed, indexGas, indexIons, NP, NG, NS, SIZE)
            % Update temp items

            % Get species to be removed
            FLAG_REMOVE = N(index, 1) / NP < exp(-SIZE);

            % Check if there are species to be removed
            if ~sum(FLAG_REMOVE)
                return
            end
            
            % Set to zero the moles of the species to be removed
            N(index(FLAG_REMOVE), 1) = 0;

            % Get the index of the species to be removed from index list
            indexRemove = find(FLAG_REMOVE)';
            
            % Update index list
            for i = indexRemove
                indexGas(indexGas == index(i)) = [];
                indexCondensed(indexCondensed == index(i)) = [];
                indexIons(indexIons == index(i)) = [];
            end

            % Update index list and lengths
            index = [indexGas, indexCondensed];
            NG = length(indexGas);
            NS = length(index);
        end

        function delta = relaxFactor(NP, ni, eta, Delta_ln_NP, NG)
            % Compute relaxation factor
            FLAG = eta(1:NG) > 0;
            FLAG_MINOR = ni(1:NG) / NP <= 1e-8 & FLAG;
            delta1 = 2./max(5*abs(Delta_ln_NP), abs(eta(FLAG)));
            delta2 = min(abs((-log(ni(FLAG_MINOR)/NP) - 9.2103404) ./ (eta(FLAG_MINOR) - Delta_ln_NP)));
            delta = min([1; delta1; delta2]);
        end

        function point = getPoint(x_vector, f_vector)
            % Get point using the regula falsi method
            %
            % Args:
            %     x_vector (float): Guess temperature [K]
            %     f_vector (struct): evaluated functions [J] (HP, EV) or [J/K] (SP, SV)
            % Returns:
            %     point (float): Point of the function [K]
        
            point = (f_vector(2) * x_vector(1) - f_vector(1) * x_vector(2)) / (f_vector(2) - f_vector(1));
        end

        function point = getPointAitken(x0, g_vector)
            % Get fixed point of a function based on the chemical transformation using the Aitken acceleration method
            %
            % Args:
            %     x0 (float): Guess temperature [K]
            %     g_vector (struct): Fixed points of the function [J] (HP, EV) or [J/K] (SP, SV)
            %
            % Returns:
            %     point (float): Point of the function [K]
            
            point = x0 - (g_vector(1) - x0)^2 / (g_vector(2) - 2*g_vector(1) + x0);   
        end

    end

    methods (Access = private)
        
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumGibbs(obj, system, p, T, mix, guess_moles)
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, v, T, mix, guess_moles)
        [N, STOP_ions, FLAG_ION] = equilibriumCheckIons(obj, N, A0, ind_E, indexGas, indexIons)

        function [x, STOP, guess_moles] = newton(obj, mix1, mix2, field, x0, guess_moles)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the second-order Newton-Raphson method
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     pP (float): Pressure [bar]
            %     field (str): Fieldname in Problem Description (PD)
            %     x0 (float): Guess temperature [K]
            %     guess_moles (float): Guess moles final mixture
            %
            % Returns:
            %     Tuple containing
            %
            %     * x (float): Temperature at equilibrium [K]
            %     * STOP (float): Relative error [-] 
            %     * guess_moles (struct): Guess moles final mixture
        
            if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                x = x0;
                STOP = 0;
                return
            end
            
            % Initialization
            it = 0; STOP = 1.0;
        
            % Loop
            while STOP > obj.tol0 && it < obj.itMax
                % Update iteration number
                it = it + 1;
                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                [f0, fprime0, frel, guess_moles] = obj.getRatioNewton(mix1, mix2, field, x0, guess_moles);
                % Compute solution
                x = abs(x0 - f0 / fprime0);
                % Compute stop criteria
                STOP = max(abs((x - x0) / x), frel);
                % Update solution
                x0 = x;
                % Debug       
                % aux_x(it) = x;
                % aux_STOP(it) = STOP;
            end
        
            % debug_plot_error(it, aux_STOP);
            if STOP > obj.tol0
                fprintf('\n***********************************************************\n')
                fprintf('Newton method not converged\nCalling Newton-Steffensen root finding algorithm\n')
                x0 = obj.regulaGuess(mix1, mix2, field);
                [x, STOP] = obj.nsteff(mix1, mix2, field, x0, []);
            else
                obj.printError(it, x, STOP);
            end
            
        end
        
        function [f, fprime, frel, guess_moles] = getRatioNewton(obj, mix1, mix2, field, x, guess_moles)
                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                
                try
                    obj.equilibrate_T(mix1, mix2, x, guess_moles);
                catch
                    obj.equilibrate_T(mix1, mix2, x);
                end
            
                % Calculate residual of f = 0
                f = mix2.(field) - mix1.(field);
                % Calculate partial derivative of f with temperature
                fprime = obj.getPartialDerivative(mix2);
                % Get relative value of the residual
                frel = abs(f / mix2.(field));
                % Update guess moles
                guess_moles = mix2.N * mix2.Xi;
        end

        function value = getPartialDerivative(obj, mix)
            % Get value of the partial derivative for the set problem type [J/K] (HP, EV) or [J/K^2] (SP, SV)
            %
            % Args:
            %     obj (EquilibriumSolver):
            %     mix (Mixture):
            %
            % Returns:
            %     value (float): Value of the partial derivative for the set problem type [J/K] (HP, EV) or [J/K^2] (SP, SV)
        
            if strcmpi(obj.problemType, 'HP')
                value = mix.cp;
            elseif strcmpi(obj.problemType, 'EV')
                value = mix.cv;
            elseif strcmpi(obj.problemType, 'SP')
                value = mix.cp / mix.T;
            elseif strcmpi(obj.problemType, 'SV')
                value = mix.cv / mix.T;
            end
        
        end

        function [x, STOP, guess_moles] = nsteff(obj, mix1, mix2, field, x0, guess_moles)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the third-order Newton-Steffensen method
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     field (str): Fieldname in Problem Description (PD)
            %     x0 (float): Guess temperature [K]
            %     guess_moles (float): Guess moles final mixture
            %
            % Returns:
            %     Tuple containing
            %
            %     * x (float): Temperature at equilibrium [K]
            %     * STOP (float): Relative error [-] 
            %     * guess_moles (struct): Guess moles final mixture
        
            if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                x = x0;
                STOP = 0;
                return
            end
        
            it = 0; STOP = 1.0;
            while STOP > obj.tol0 && it < obj.itMax
                % Update iteration number
                it = it + 1;
                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                [f0, fprime0, frel, guess_moles] = getRatioNewton(obj, mix1, mix2, field, x0, guess_moles);
                % Compute pseudo-solution
                x = abs(x0 - f0 / fprime0);
                % Re-estimation of first derivative
                f0_2 = getRatioNewton(obj, mix1, mix2, field, x, guess_moles);
                % Compute solution
                x = abs(x0 - f0^2 / (fprime0 * (f0 - f0_2)));
                % Compute stop criteria
                STOP = max(abs((x - x0) / x), frel);
                % Update solution
                x0 = x;
                % Debug       
                % aux_x(it) = x;
                % aux_STOP(it) = STOP;
            end

            % debug_plot_error(it, aux_STOP);

            if STOP > obj.tol0
                fprintf('\n***********************************************************\n')
                fprintf('Newton method not converged\nCalling Steffensen-Aitken root finding algorithm\n')
                x0 = obj.regulaGuess(mix1, mix2, field);
                [x, STOP] = obj.steff(mix1, mix2, field, x0, []);
            else
                obj.printError(it, x, STOP);
            end
    
        end

        function [x, STOP, guess_moles] = steff(obj, mix1, mix2, field, x0, guess_moles)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the Steffenson-Aitken method
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     field (str): Fieldname in Problem Description (PD)
            %     x0 (float): Guess temperature [K]
            %     guess_moles (float): Guess moles final mixture
            %
            % Returns:
            %     Tuple containing
            %
            %     * x (float): Temperature at equilibrium [K]
            %     * STOP (float): Relative error [-] 
            %     * guess_moles (float): Guess moles final mixture
        
            if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                x = x0;
                STOP = 0;
                return
            end
        
            it = 0; STOP = 1.0;
            
            while STOP > obj.tol0 && it < obj.itMax
                it = it + 1;
                [g, g_rel]= getGpoint(obj, mix1, mix2, field, x0, guess_moles);
                fx = abs(g - x0);
                g_aux  = getGpoint(obj, mix1, mix2, field, fx, guess_moles);
                fx2 = abs(g_aux - fx);
                if abs(fx2 - 2*fx + x0) > obj.tol0
                    x = obj.getPointAitken(x0, [fx, fx2]);
                else
                    x = fx;
                end
        
                STOP = max(abs((x - fx) / x), abs(g_rel));
                x0 = x;
            end
        
            obj.printError(it, x, STOP);
        end

        function x0 = regulaGuess(obj, mix1, mix2, field)
            % Find a estimate of the temperature for the set chemical equilibrium
            % transformation using the regula falsi method
            %
            % Args:
            %     obj (EquilibriumSolver): 
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     field (str): Fieldname in Problem Description (PD)
            %
            % Returns:
            %     x0 (float): Guess temperature [K]
            
            % Initialization
            guess_moles = [];

            % Define branch
            x_l = obj.root_T0_l;
            x_r = obj.root_T0_r;
            
            % Compute f(x) = f2(x) - f1 at the branch limits
            g_l = obj.getGpoint(mix1, mix2, field, x_l, guess_moles);
            g_r = obj.getGpoint(mix1, mix2, field, x_r, guess_moles);
            
            % Update estimate based on the region
            if g_l * g_r < 0
                x0 = obj.getPoint([x_l, x_r], [g_l, g_r]);
            elseif g_l * g_r > 0 && abs(g_l) < abs(g_r)
                x0 = x_l - 50;
            elseif g_l * g_r > 0 || (isnan(g_l) && isnan(g_r))
                x0 = x_r + 50;
            elseif isnan(g_l) && ~isnan(g_r)
                x0 = x_r - 100;
            elseif ~isnan(g_l) && isnan(g_r)
                x0 = x_l + 100;
            else
                x0 = obj.root_T0;
            end
        
        end

        function [gpoint, gpoint_relative] = getGpoint(obj, mix1, mix2, field, x0, guess_moles)
            % Get fixed point of a function based on the chemical transformation
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     field (str):    Fieldname in Problem Description (PD)
            %     x0 (float):     Guess temperature [K]
            %
            % Returns:
            %     Tuple containing
            %
            %     * gpoint (float): Fixed point of the function [J] (HP, EV) or [J/K] (SP, SV)
            %     * gpoint_relative (float): Fixed relative point of the function [J] (HP, EV) or [J/K] (SP, SV)
            
            try
                % Compute TP problem
                obj.equilibrate_T(mix1, mix2, x0, guess_moles);
                % Compute f(x) = f2(x) - f1 = 0
                gpoint = (mix2.(field) - mix1.(field));
                % Compute f(x) / f2(x)
                gpoint_relative = gpoint / (mix2.(field));
        
                if strcmpi(field, 's')
                    gpoint = gpoint * 1e3;
                end
        
            catch
                gpoint = NaN;
                gpoint_relative = NaN;
            end

        end

        function printError(obj, it, T, STOP)
            % Print error of the method if the number of iterations is greater than maximum iterations allowed
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     it (float):    Number of iterations executed in the method
            %     T (float):     Temperature [K]
            %     STOP (float):  Relative error [-]
        
            if it < obj.itMax
                return
            end

            fprintf('***********************************************************\n')
            fprintf('Root algorithm not converged \n')
            fprintf('   Error       =  %8.2f [%%]  \n', STOP*100)
            fprintf('   Temperature =  %8.2f [K]  \n', T)
            fprintf('   Iterations  =  %8.d [it] \n', it)
            fprintf('***********************************************************\n')
        end
    
    end

end