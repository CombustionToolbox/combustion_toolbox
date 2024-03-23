classdef EquilibriumSolver
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
        root_method = @newton      % Root finding method [newton (2nd order), steff (2nd order), or nsteff (3rd order)]
        root_T0_l = 1000           % First temperature guess [K] left branch - root finding method
        root_T0_r = 3000           % First temperature guess [K] right branch - root finding method
        root_T0   = 3000           % Temperature guess [K] if it's outside previous range - root finding method
        % * Flags
        FLAG_EXTRAPOLATE = true;   % Flag indicating linear extrapolation of the polynomials fits
        FLAG_FAST = true;          % Flag indicating use guess composition of the previous computation
        FLAG_TCHEM_FROZEN = false; % Flag to consider a thermochemically frozen gas (calorically perfect gas)
        FLAG_FROZEN = false;       % Flag to consider a calorically imperfect gas with frozen chemistry
        FLAG_EOS = false;          % Flag to use non-ideal Equation of States (EoS)
    end

    methods

        function obj = EquilibriumSolver(varargin)
            % Constructor
            defaultProblemType = 'TP';

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'TP', 'TV', 'HP', 'EV', 'SP', 'SV'})));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
        end

        function mix = solve(obj, mix)
            % Obtain chemical equilibrium composition and thermodynamic properties
            mix = obj.equilibrate(mix);

        end

        function mix = equilibrate(obj, mix)
            % Obtain properties at equilibrium for the given thermochemical transformation
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            %     mix (Mixture): Mixture considering a thermochemical frozen gas
            %
            % Returns:
            %     mix (struct): Mixture at chemical equilibrium for the given thermochemical transformation
            %
            % Example:
            %     * obj = EquilibriumSolver();
            %     * mix = obj.equilibrate('TP', mix)

            if obj.FLAG_TCHEM_FROZEN
                mix.errorProblem = 0;
                return
            end

            % get attribute xx of the specified transformations
            attr_name = get_attr_name(obj);

            % compute initial guess
            [guess, guess_moles] = get_guess(obj, mix, attr_name);

            % If the problem type is SV, the product's volume is based on the given v_P/v_R ratio
            mix = set_volume_SV(obj, mix);

            % root finding: find the value x that satisfies f(x) = mix2.xx(x) - mix1.xx = 0
            [T, STOP, guess_moles] = root_finding(obj, mix, attr_name, guess, guess_moles);

            % compute properties
            mix = equilibrate_T(obj, mix, T, guess_moles);

            % check convergence in case the problemType is TP (defined Temperature and Pressure)
            print_convergence(mix.errorMoles, obj.tolGibbs, mix.errorMolesIons, obj.tolMultiplierIons, obj.problemType)

            % save error from root finding algorithm
            mix.errorProblem = STOP;

            % save equivalence ratio
            % mix2.equivalenceRatio = get_phi(mix1);

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

            function [guess, guess_moles] = get_guess(obj, mix, attr_name)
                % Get initial estimates for temperature and molar composition

                % Initialization
                if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                    % guess = get_transformation(obj, 'TP');
                    guess = mix.T;
                    if mix.FLAG_REACTION
                        guess_moles = mix.Xi * mix.N;
                    else
                        guess_moles = [];
                    end
                    
                elseif ~isempty(mix)
                    guess = mix.T;
                    guess_moles = mix.Xi * mix.N;
                else
                    guess = regula_guess(obj, mix1, pP, attr_name);
                    guess_moles = [];
                end

            end

            function [x, STOP, guess_moles] = root_finding(obj, mix, attr_name, x0, guess_moles)
                % Calculate the temperature value that satisfied the problem conditions
                % using the @root_method
                [x, STOP, guess_moles] = obj.root_method(obj, mix, attr_name, x0, guess_moles);
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

        function mix = equilibrate_T(obj, mix, TP, varargin)
            % Obtain equilibrium properties and composition for the given temperature [K] and pressure [bar]
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %     mix1 (struct): Properties of the initial mixture
            %     pP (float): Pressure [bar]
            %     TP (float): Temperature [K]
            %
            % Optional Args:
            %     guess_moles (float): Mixture composition [mol] of a previous computation
            %
            % Returns:
            %     mix2 (struct): Properties of the final mixture
            %
            % Example:
            %     mix2 = equilibrate(self, self.PS.strR{1}, 1.01325, 3000)
            
            % Check if calculations are for a thermochemical frozen gas (calorically perfect)
            if obj.FLAG_TCHEM_FROZEN
                return
            end
        
            % Check if calculations are for a calorically imperfect gas with frozen chemistry
            if obj.FLAG_FROZEN
                mix = equilibrate_T_frozen(obj, mix1, pP, TP);
                return
            end
            
            % Definitions
            N_mix0 = moles(mix); % Get moles of inert species

            % Unpack
            guess_moles = unpack(varargin);
            % Check flag
            if ~obj.FLAG_FAST, guess_moles = []; end
            % Set List of Species to List of Products
            system = set_LS_original(obj, mix.chemicalSystem, TP);
            % Compute number of moles
            [N, dNi_T, dN_T, dNi_p, dN_p, indexListProducts, STOP, STOP_ions, h0] = select_equilibrium(obj, system, TP, mix, guess_moles);
            % Get index of species
            index = find_ind(system.listSpecies, system.listSpecies(indexListProducts));
            % Reshape composition matrix N, and partial composition partial derivatives 
            N = reshape_vector(system, index, indexListProducts, N);
            dNi_T = reshape_vector(system, index, indexListProducts, dNi_T);
            dNi_p = reshape_vector(system, index, indexListProducts, dNi_p);
            N(system.indexFrozen) = N_mix0(system.indexFrozen);
            % Compute property matrix of the species at chemical equilibrium
            % NOTE: If the ind variable is removed from the inputs, the set_species 
            % routine will completely fill the properties matrix
            % M0 = set_species(obj, system.listSpecies, N(:, 1), TP, [index, system.indexFrozen]);
            
            mix.chemicalSystem = system.clean;
            mix.T = TP;
            mix.dNi_T = dNi_T; mix.dN_T = dN_T;
            mix.dNi_p = dNi_p; mix.dN_p = dN_p;
            mix.FLAG_REACTION = true;
            mix.errorMoles = STOP;
            mix.errorMolesIons = STOP_ions;
            mix = mix.set_fast(system.listSpecies, N(:, 1)', [index, system.indexFrozen], h0);

            % Compute properties of final mixture
            % mix2 = set_properties(obj, mix1, M0, pP, TP, STOP, STOP_ions);

            % SUB-PASS FUNCTIONS
            function guess_moles = unpack(value)
                % Unpack inputs
                if isempty(value)
                    guess_moles = [];
                else
                    guess_moles = value{1};
                end
            
            end
            
            function system = set_LS_original(obj, system, TP)
                % Set List of Species to List of Products
                
                % Remove ionized species if TP is below T_ions
                if any(system.indexIons) && TP < obj.temperatureIons
                    system.indexListSpeciesOriginal(system.indexIons) = [];
                    system.indexElements = [];
                end
                
                % Initialization
                system.indexGas = []; system.indexCondensed = [];
                system.indexCryogenic = []; system.indexIons = [];
                % Set list of species for calculations
                system.listSpecies = system.listSpecies(system.indexListSpeciesOriginal);
                % Establish cataloged list of species according to the state of the phase
                system = system.list_phase_species();
                % Update stoichiometric matrix
                system.stoichiometricMatrix = system.stoichiometricMatrix(system.indexListSpeciesOriginal, :);
                % Update property matrix
                system.propertiesMatrix = system.propertiesMatrix(system.indexListSpeciesOriginal, :);
                % Update compostion matrix
                system.molesPhaseMatrix = system.molesPhaseMatrix(system.indexListSpeciesOriginal, :);
            end
            
            function pP = set_pressure(self, mix1, TP, N)
                % Compute pressure of product mixture
                pP = self.PD.EOS.pressure(self, N, TP, mix1.v, system.LS, mix1.Xi) * 1e-5;
            end
            
            function [N, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions, h0] = select_equilibrium(obj, system, TP, mix, guess_moles)
                % Select equilibrium: TP: Gibbs; TV: Helmholtz
                if strfind(obj.problemType, 'P') == 2
                    [N, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions, h0] = equilibrium_gibbs(obj, system, mix.p, TP, mix, guess_moles);
                else
                    [N, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions] = equilibrium_helmholtz(obj, system, mix.v, TP, mix, guess_moles);
                end

                h0 = h0 * 1e-3; % [kJ/mol]
            end
            
            function mix2 = set_properties(self, mix1, M0, pP, TP, STOP, STOP_ions)
                % Compute properties of final mixture
                if strfind(self.PD.ProblemType, 'P') == 2
                    mix2 = compute_properties(self, M0, pP, TP);
                else
                    NP = sum(M0(:, self.C.M0.ind_ni) .* (1 - M0(:, self.C.M0.ind_phase)));
                    pP = set_pressure(self, mix1, TP, NP);
                    mix2 = compute_properties(self, M0, pP, TP);
                end
            
                mix2.error_moles = STOP;
                mix2.error_moles_ions = STOP_ions;
            end
            
            function vector = reshape_vector(system, index, ind_modified, vector_modified)
                % Reshape vector containing all the species
                vector = system.molesPhaseMatrix(:, 1);
                vector(index, 1) = vector_modified(ind_modified, 1);
            end
            
            function mix2 = equilibrate_T_frozen(self, mix1, pP, TP)
                % Obtain equilibrium properties and composition for the given
                % temperature [K] and pressure [bar] assuming a calorically imperfect
                % gas with frozen chemistry
                
                % Initialization
                N = self.C.N0.value; % Composition matrix [n_i, FLAG_CONDENSED_i]
            
                % Get data from mix1
                STOP = mix1.error_moles;
                STOP_ions = mix1.error_moles_ions;
            
                % Set all species as frozen
                self = set_react_index(self, system.LS(self.Misc.index_LS_original));
                index = system.indexFrozen;
            
                % Add moles of frozen species to the moles vector N
                N_mix1 = moles(mix1);
                N(system.indexFrozen) = N_mix1(system.indexFrozen);
            
                % Compute property matrix of the species at chemical equilibrium
                % NOTE: If the ind variable is removed from the inputs, the set_species 
                % routine will completely fill the properties matrix
                M0 = set_species(self, system.LS, N(:, 1), TP, index);
            
                % Compute properties of final mixture
                mix2 = set_properties(self, mix1, M0, pP, TP, STOP, STOP_ions);
            
                % Set thermodynamic derivative to their frozen values
                mix2.dVdT_p = 1;
                mix2.dVdp_T = -1;
            end

        end

    end

    methods (Access = private)
        
        function [x, STOP, guess_moles] = newton(obj, mix, field, x0, guess_moles)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the second-order Newton-Raphson method
            %
            % Args:
            %     self (struct): Data of the mixture, conditions, and databases
            %     mix1 (struct): Properties of the initial mixture
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
            
            % Definitions
            pP = mix.p;

            % Initialization
            it = 0; STOP = 1.0;
        
            % Loop
            while STOP > obj.tol0 && it < obj.itMax
                % Update iteration number
                it = it + 1;
                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                [f0, fprime0, frel, guess_moles] = get_ratio_newton(obj, mix, field, x0, guess_moles);
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
            if STOP > obj.TN.tol0
                fprintf('\n***********************************************************\n')
                fprintf('Newton method not converged\nCalling Newton-Steffensen root finding algorithm\n')
                x0 = regula_guess(obj, mix, pP, field);
                [x, STOP] = nsteff(obj, mix, pP, field, x0, []);
            else
                print_error_root(it, obj.TN.itMax, x, STOP);
            end

            % SUB-PASS FUNCTIONS
            function [f, fprime, frel, guess_moles] = get_ratio_newton(obj, mix, field, x, guess_moles)
                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                try
                    mix = equilibrate_T(obj, mix, x, guess_moles);
                catch
                    mix = equilibrate_T(obj, mix, x);
                end
            
                % Calculate residual of f = 0
                f = mix.(field) - mix.(field);
                % Calculate partial derivative of f with temperature
                fprime = get_partial_derivative(obj, mix);
                % Get relative value of the residual
                frel = abs(f / mix.(field));
                % Update guess moles
                guess_moles = mix.N * mix.Xi;
            end
            
        end
    
    end

end