classdef EquilibriumSolver < handle
    % EquilibriumSolver: A class for solving chemical equilibrium problems
    
    properties
        problemType                % Problem type [TP, TV, HP, EV, SP, SV]
        % * Chemical equilibrium TP, TV (CT-EQUIL module)
        tolGibbs = 1e-5            % Tolerance of the Gibbs/Helmholtz minimization method
        tolE = 1e-6                % Tolerance of the mass balance
        tolMoles = 1e-14           % Tolerance of the composition of the mixture                       
        tolMolesGuess = 1e-6       % Tolerance of the molar composition of the mixture (guess)
        tolMultiplierIons = 1e-4   % Tolerance of the dimensionless Lagrangian multiplier - ions
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
        FLAG_EXTRAPOLATE = true;   % Flag indicating linear extrapolation of the polynomials fits
        FLAG_FAST = true;          % Flag indicating use guess composition of the previous computation
        FLAG_TCHEM_FROZEN = false; % Flag to consider a thermochemically frozen gas (calorically perfect gas)
        FLAG_FROZEN = false;       % Flag to consider a calorically imperfect gas with frozen chemistry
        FLAG_EOS = false;          % Flag to use non-ideal Equation of States (EoS)
        % Miscellaneous
        FLAG_RESULTS = true        % Flag to print results
    end

    methods

        function obj = EquilibriumSolver(varargin)
            % Constructor
            defaultProblemType = 'TP';

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'TP', 'TV', 'HP', 'EV', 'SP', 'SV'})));
            addOptional(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
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

            % Calculations
            obj.solve(mixArray(n));
            
            for i = n-1:-1:1
                obj.solve(mixArray(i), mixArray(i + 1));
            end

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
            mix2 = set_volume_SV(obj, mix2);
            % Root finding: find the value x that satisfies f(x) = mix2.xx(x) - mix1.xx = 0
            [T, STOP, guess_moles] = rootFinding(obj, mix1, mix2, attr_name, guess, guess_moles);
            % Compute properties
            obj.equilibrate_T(mix2, T, guess_moles);
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

                mix.v = mix.v * obj.PD.vP_vR.value;
            end

        end

        function mix = equilibrate_T(obj, mix, T, varargin)
            % Obtain equilibrium properties and composition for the given temperature [K] and pressure [bar]
            %
            % Args:
            %     obj (EquilibriumSolver): 
            %     mix (Mixture): Properties of the initial mixture
            %     T (float): Temperature [K]
            %
            % Optional Args:
            %     guess_moles (float): Mixture composition [mol] of a previous computation
            %
            % Returns:
            %     mix (Mixture): Properties of the final mixture
            %
            % Example:
            %     mix = equilibrate(EquilibriumSolver, mix, 3000)
            
            % Import packages
            import combustiontoolbox.utils.findIndex
            
            % Check if calculations are for a thermochemical frozen gas (calorically perfect)
            if obj.FLAG_TCHEM_FROZEN
                % TO BE IMPLEMENTED
                return
            end
        
            % Check if calculations are for a calorically imperfect gas with frozen chemistry
            if obj.FLAG_FROZEN
                % Computed by default when defining the mixture (Mixture)
                return
            end
            
            % Definitions
            N_mix0 = moles(mix); % Get moles of inert species
            system = mix.chemicalSystem;
            systemProducts = mix.chemicalSystemProducts;
            % Unpack
            guess_moles = unpack(varargin);
            % Check flag
            if ~obj.FLAG_FAST, guess_moles = []; end
            % Compute number of moles
            [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = selectEquilibrium(obj, systemProducts, T, mix, guess_moles);
            % Compute property matrix of the species at chemical equilibrium
            mix = setMixture(mix);

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
                vMolar = vSpecific2vMolar(mix, mix.vSpecific, moles, index);
                pP = mix.equationOfState.getPressure(T, vMolar, system.listSpecies, mix.Xi) * 1e-5;
            end
            
            function [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = selectEquilibrium(obj, system, T, mix, guess_moles)
                % Select equilibrium: TP: Gibbs; TV: Helmholtz
                
                if strfind(obj.problemType, 'P') == 2
                    [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumGibbs(obj, system, mix.p, T, mix, guess_moles);
                else
                    [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, mix.v, T, mix, guess_moles);
                end

                h0 = h0 * 1e-3; % [kJ/mol]
            end
            
            function mix = setMixture(mix)
                % Compute properties of final mixture
                
                % Reshape composition matrix N, and partial composition partial derivatives 
                N = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, N);
                dNi_T = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, dNi_T);
                dNi_p = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, dNi_p);
                h0 = reshapeVector(system, system.indexProducts, systemProducts.indexSpecies, h0);
                % h0(system.indexFrozen) = set_h0(system.listSpecies(system.indexFrozen), T, system.species) * 1e-3; 
                N(system.indexFrozen) = N_mix0(system.indexFrozen);

                % Assign values
                mix.T = T;
                mix.dNi_T = dNi_T; mix.dN_T = dN_T;
                mix.dNi_p = dNi_p; mix.dN_p = dN_p;
                mix.FLAG_REACTION = true;
                mix.errorMoles = STOP;
                mix.errorMolesIons = STOP_ions;

                % Clean system
                system.clean;
                
                % Check if problemType is at constant volume
                if strfind(obj.problemType, 'V') == 2
                    mix.p = computePressure(mix, T, N, system.indexProducts);
                end
                
                % Compute properties of final mixture
                mix.set_fast(system.listSpecies, N', [indexProducts, system.indexFrozen], h0);
            end
            
            function vector = reshapeVector(system, index, ind_modified, vector_modified)
                % Reshape vector containing all the species
                vector = system.molesPhaseMatrix(:, 1);
                vector(index, 1) = vector_modified(ind_modified, 1);
            end
            
        end

    end
    
    methods (Access = private, Static)
        [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N0, A0, NE, indexGas, indexCondensed, indexElements, H0RT);

        function point = getPoint(x_vector, f_vector)
            % Get point using the regula falsi method
            %
            % Args:
            %     x_vector (float): Guess temperature [K]
            %     f_vector (struct): evaluated functions [kJ] (HP, EV) or [kJ/K] (SP, SV)
            % Returns:
            %     point (float): Point of the function [K]
        
            point = (f_vector(2) * x_vector(1) - f_vector(1) * x_vector(2)) / (f_vector(2) - f_vector(1));
        end

    end

    methods (Access = private)
        
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumGibbs(obj, system, p, T, mix, guess_moles)
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, vP, T, mix, guess_moles)

        function [x, STOP, guess_moles] = newton(obj, mix1, mix2, field, x0, guess_moles)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the second-order Newton-Raphson method
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %     mix1 (struct): Properties of the initial mixture
            %     mix2 (struct): Properties of the final mixture
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
                    obj.equilibrate_T(mix2, x, guess_moles);
                catch
                    obj.equilibrate_T(mix2, x);
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
            % Get value of the partial derivative for the set problem type [kJ/K] (HP, EV) or [kJ/K^2] (SP, SV)
            %
            % Args:
            %     obj (EquilibriumSolver):
            %     mix (Mixture):
            %
            % Returns:
            %     value (float): Value of the partial derivative for the set problem type [kJ/K] (HP, EV) or [kJ/K^2] (SP, SV)
        
            if strcmpi(obj.problemType, 'HP')
                value = mix.cp;
            elseif strcmpi(obj.problemType, 'EV')
                value = mix.cv;
            elseif strcmpi(obj.problemType, 'SP')
                value = mix.cp / mix.T;
            elseif strcmpi(obj.problemType, 'SV')
                value = mix.cv / mix.T;
            end
        
            value = value * 1e-3; % [kJ/K] (HP, EV) or [kJ/K^2] (SP, SV)
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
            %     self (struct):  Data of the mixture, conditions, and databases
            %     mix1 (struct):  Properties of the initial mixture
            %     mix2 (struct):  Properties of the final mixture
            %     field (str):    Fieldname in Problem Description (PD)
            %     x0 (float):     Guess temperature [K]
            %
            % Returns:
            %     Tuple containing
            %
            %     * gpoint (float): Fixed point of the function [kJ] (HP, EV) or [kJ/K] (SP, SV)
            %     * gpoint_relative (float): Fixed relative point of the function [kJ] (HP, EV) or [kJ/K] (SP, SV)
            
            try
                % Compute TP problem
                obj.equilibrate_T(mix2, x0, guess_moles);
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
            %     it (float):    Number of iterations executed in the method
            %     itMax (float): Maximum nNumber of iterations allowed in the method
            %     T (float):     Temperature [K]
            %     STOP (float):  Relative error [-]
        
            if it == obj.itMax
                fprintf('***********************************************************\n')
                fprintf('Root algorithm not converged \n')
                fprintf('   Error       =  %8.2f [%%]  \n', STOP*100)
                fprintf('   Temperature =  %8.2f [K]  \n', T)
                fprintf('   Iterations  =  %8.d [it] \n', it)
                fprintf('***********************************************************\n')
            end

        end
    
    end

end