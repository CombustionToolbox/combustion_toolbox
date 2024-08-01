classdef EquilibriumSolver < handle
    % The :mat:func:`EquilibriumSolver` class is used to compute the composition at the equilibrium
    % of multi-component gas mixtures that undergo canonical thermochemical transformations from
    % an initial state (reactants), defined by its initial composition, temperature, and pressure,
    % to a final state (products), defined by a set of chemical species (in gaseous---included
    % ions---or pure condensed phase). 
    % 
    % The :mat:func:`EquilibriumSolver` object can be initialized as follows: ::
    %
    %       solver = EquilibriumSolver('problemType', problemType, ...)
    %
    % Here ``problemType`` represents the acronym of the problem to be solved (see below).
    % Additional optional parameters can be provided to customize the solver's behavior.
    %
    % Problem types:
    %     * `TP`: Equilibrium composition at defined temperature and pressure
    %     * `TV`: Equilibrium composition at defined temperature and specific volume
    %     * `HP`: Equilibrium composition at defined enthalpy and pressure
    %     * `EV`: Equilibrium composition at defined internal energy and specific volume
    %     * `SP`: Equilibrium composition at defined entropy and pressure
    %     * `SV`: Equilibrium composition at defined entropy and specific volume
    %
    % See also: :mat:func:`Mixture`, :mat:func:`solve`, :mat:func:`solveArray`, :mat:func:`report`
    
    properties
        problemType                % Problem type [TP, TV, HP, EV, SP, SV]
        tolGibbs = 1e-6            % Tolerance of the Gibbs/Helmholtz minimization method
        tolE = 1e-6                % Tolerance of the mass balance
        tolMoles = 1e-14           % Tolerance of the composition of the mixture                       
        tolMolesGuess = 1e-6       % Tolerance of the molar composition of the mixture (guess)
        tolMultiplierIons = 1e-4   % Tolerance of the dimensionless Lagrangian multiplier - ions
        tolTau = 1e-25             % Tolerance of the slack variables for condensed species
        itMaxGibbs = 70            % Max number of iterations - Gibbs/Helmholtz minimization method
        itMaxIons = 30             % Max number of iterations - charge balance (ions)
        temperatureIons = 0        % Minimum temperature [K] to consider ionized species
        tol0 = 1e-3                % Tolerance of the root finding algorithm
        itMax = 30                 % Max number of iterations - root finding method
        rootMethod = @newton       % Root finding method [newton (2nd order), steff (2nd order), or nsteff (3rd order)]
        root_T0_l = 1000           % First temperature guess [K] left branch - root finding method
        root_T0_r = 3000           % First temperature guess [K] right branch - root finding method
        root_T0   = 3000           % Temperature guess [K] if it's outside previous range - root finding method
        FLAG_EXTRAPOLATE = true    % Flag indicating linear extrapolation of the polynomials fits
        FLAG_FAST = true           % Flag indicating use guess composition of the previous computation
        FLAG_TCHEM_FROZEN = false  % Flag to consider a thermochemically frozen gas (calorically perfect gas)
        FLAG_FROZEN = false        % Flag to consider a calorically imperfect gas with frozen chemistry
        FLAG_EOS = false           % Flag to use non-ideal Equation of States (EoS)
        FLAG_RESULTS = true        % Flag to print results
        FLAG_TIME = true           % Flag to print elapsed time
        FLAG_REPORT = false        % Flag to postprocess all the results with predefined plots
        time                       % Elapsed time [s]
        plotConfig                 % PlotConfig object
    end

    methods
        
        mix2 = equilibrate(obj, mix2, varargin)
        mix2 = equilibrateT(obj, mix1, mix2, T, varargin)

        function obj = EquilibriumSolver(varargin)
            % Constructor
            defaultProblemType = 'TP';
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();

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
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
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
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.plotConfig = p.Results.plotConfig;
        end

        function mix = solve(obj, mix, varargin)
            % Obtain chemical equilibrium composition and thermodynamic properties
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix (Mixture): Mixture object containing initial state
            %
            % Optional Args:
            %     * mixGuess (Mixture): Mixture object from previous calculation
            %
            % Returns:
            %     mix (Mixture): Mixture object updated with equilibrium composition and properties
            %
            % Examples:
            %     * mix = solve(EquilibriumSolver(), mix);
            %     * mix = solve(EquilibriumSolver(), mix, mixGuess);

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
            % Obtain chemical equilibrium composition and thermodynamic properties for an array of mixture values
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mixArray (array of Mixture): Array of Mixture objects containing initial states
            %
            % Returns:
            %     mixArray (array of Mixture): Array of Mixture objects updated with equilibrium compositions and properties
            %
            % Example:
            %     * mixArray = solveArray(EquilibriumSolver(), mixArray);
            
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

            % Postprocess all the results with predefined plots
            if obj.FLAG_REPORT
                report(obj, mixArray);
            end

        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            
            if ~obj.FLAG_TIME
                return
            end

            fprintf('\nElapsed time is %.5f seconds\n', obj.time);
        end

        function ax2 = plot(obj, mixArray, varargin)
            % Plot results
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mixArray (Mixture): Array of Mixture objects
            %
            % Optional Args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Returns:
            %     ax2 (Axes): Axes object with properties of the mixtures
            %
            % Examples:
            %     * plot(EquilibriumSolver(), mixArray);
            %     * plot(EquilibriumSolver(), mixArray, mixArray2);
            
            % Import packages
            import combustiontoolbox.utils.display.plotComposition
            import combustiontoolbox.utils.display.plotProperties
            
            % Definitions
            additionalMixtures = nargin - 2;
            numPlotProperties = obj.plotConfig.numPlotProperties;

            % Check if is a scalar value
            if isscalar(mixArray)
                return
            end

            % Get labels
            if additionalMixtures

                if mixArray(1).chemicalSystem.FLAG_COMPLETE && ~varargin{1}(1).chemicalSystem.FLAG_COMPLETE
                    labels = {'Complete', 'Incomplete'};
                else
                    labels = arrayfun(@(x) sprintf('Mixture %d', x), 1:(additionalMixtures + 1), 'UniformOutput', false);
                end

            end
            
            % Plot molar fractions - mixArray
            ax1 = plotComposition(mixArray(1), mixArray, mixArray(1).rangeName, 'Xi', 'mintol', obj.plotConfig.mintolDisplay);
        
            % Plot properties - mixArray
            ax2 = plotProperties(repmat({mixArray(1).rangeName}, 1, numPlotProperties), mixArray, obj.plotConfig.plotProperties, mixArray, 'basis', obj.plotConfig.plotPropertiesBasis, 'config', obj.plotConfig);
            
            % Check if there are additional mixtures
            if ~additionalMixtures
                return
            end

            for i = 1:nargin - 2
                % Unpack input
                mixArray = varargin{i};

                % Plot molar fractions - mixArray_i
                ax1 = plotComposition(mixArray(1), mixArray, mixArray(1).rangeName, 'Xi', 'mintol', obj.plotConfig.mintolDisplay);
            
                % Plot properties - mixArray_i
                ax2 = plotProperties(repmat({mixArray(1).rangeName}, 1, numPlotProperties), mixArray, obj.plotConfig.plotProperties, mixArray, 'basis', obj.plotConfig.plotPropertiesBasis, 'config', obj.plotConfig, 'ax', ax2);
            end

            % Set legends
            legend(ax2.Children(end), labels, 'Interpreter', 'latex', 'FontSize', ax2.Children(end).FontSize);
        end

        function ax = report(obj, mixArray, varargin)
            % Postprocess all the results with predefined plots
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mixArray (Mixture): Array of Mixture objects
            %
            % Optional args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Returns:
            %     ax (Axes): Axes object with properties of the mixtures
            %
            % Examples:
            %     * ax = report(EquilibriumSolver(), mixArray);
            %     * ax = report(EquilibriumSolver(), mixArray, mixArray2);

            if nargin > 2
                ax = obj.plot(mixArray, varargin{:});
            else
                ax = obj.plot(mixArray);
            end

        end

    end
    
    methods (Access = private, Static)
        [N, NP] = equilibriumGuess(N, NP, A0, muRT, b0, index, indexGas, indexIons, NG, molesGuess)
        [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N, A0, NE, indexGas, indexCondensed, indexElements, H0RT)
        [indexCondensed, FLAG_CONDENSED, dL_dnj] = equilibriumCheckCondensed(A0, pi_i, W, indexCondensed, muRT, NC_max, FLAG_ONE, FLAG_RULE)

        function [A0, indexRemoveSpecies, ind_E, NatomE] = removeElements(NatomE, A0, ind_E, tol)
            % Find zero sum elements and remove them from stoichiometrix matrix
            %
            % Args:
            %     NatomE (float): Number of atoms of each element
            %     A0 (float): Stoichiometrix matrix
            %     ind_E (float): Index of electron element
            %     tol (float): Tolerance
            %
            % Returns:
            %     * A0 (float): Stoichiometrix matrix updated
            %     * indexRemoveSpecies (float): Index of species to be removed
            %     * ind_E (float): Index of electron element updated
            %     * NatomE (float): Number of atoms of each element updated
        
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
            
            % Check position electron "element"
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
            % Update temporal values

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
            %     x0 (float): Initial guess for temperature [K]
            %     g_vector (struct): Fixed points of the function [J] (HP, EV) or [J/K] (SP, SV)
            %
            % Returns:
            %     point (float): Point of the function [K]
            
            point = x0 - (g_vector(1) - x0)^2 / (g_vector(2) - 2*g_vector(1) + x0);   
        end

    end

    methods (Access = private)
        
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumGibbs(obj, system, p, T, mix, molesGuess)
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, v, T, mix, molesGuess)
        [N, STOP_ions, FLAG_ION] = equilibriumCheckIons(obj, N, A0, ind_E, indexGas, indexIons)

        function [x, STOP, molesGuess] = newton(obj, mix1, mix2, attributeName, x0, molesGuess)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the second-order Newton-Raphson method
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     pP (float): Pressure [bar]
            %     attributeName (char): Attribute name of the problem type
            %     x0 (float): Initial guess for temperature [K]
            %     molesGuess (float): Initial guess for moles in the final mixture
            %
            % Returns:
            %     Tuple containing
            %
            %     * x (float): Temperature at equilibrium [K]
            %     * STOP (float): Relative error [-] 
            %     * molesGuess (float): Updated guess for moles in the final mixture
        
            if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                x = x0;
                STOP = 0;
                return
            end
            
            % Initialization
            it = 0; STOP = 1.0;
        
            % Loop
            while STOP > obj.tol0 && it < obj.itMax
                % Update iteration
                it = it + 1;

                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                [f0, fprime0, frel, molesGuess] = obj.getRatioNewton(mix1, mix2, attributeName, x0, molesGuess);
                
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

            % Debug
            % debug_plot_error(it, aux_STOP);
            
            if STOP > obj.tol0
                fprintf('\n***********************************************************\n')
                fprintf('Newton method not converged\nCalling Newton-Steffensen root finding algorithm\n')
                x0 = regulaGuess(obj, mix1, mix2, attributeName);
                [x, STOP] = nsteff(obj, mix1, mix2, attributeName, x0, []);
                return
            end
            
            printError(obj, it, x, STOP);
        end
        
        function [f, fprime, frel, molesGuess] = getRatioNewton(obj, mix1, mix2, attributeName, x, molesGuess)
            % Get the residual of f, its derivative with temperature, and the
            % relative value of the residual
            
            try
                obj.equilibrateT(mix1, mix2, x, molesGuess);
            catch
                obj.equilibrateT(mix1, mix2, x);
            end
        
            % Calculate residual of f = 0
            f = mix2.(attributeName) - mix1.(attributeName);

            % Calculate partial derivative of f with temperature
            fprime = obj.getPartialDerivative(mix2);

            % Get relative value of the residual
            frel = abs(f / mix2.(attributeName));

            % Update guess moles
            molesGuess = mix2.N * mix2.Xi;
        end

        function value = getPartialDerivative(obj, mix)
            % Get value of the partial derivative for the set problem type [J/K] (HP, EV) or [J/K^2] (SP, SV)
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix (Mixture): Mixture object
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

        function [x, STOP, molesGuess] = nsteff(obj, mix1, mix2, attributeName, x0, molesGuess)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the third-order Newton-Steffensen method
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     attributeName (char): Attribute name of the problem type
            %     x0 (float): Initial guess for temperature [K]
            %     molesGuess (float): Initial guess for moles in the final mixture
            %
            % Returns:
            %     Tuple containing
            %
            %     * x (float): Temperature at equilibrium [K]
            %     * STOP (float): Relative error [-] 
            %     * molesGuess (float): Updated guess for moles in the final mixture
        
            if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                x = x0;
                STOP = 0;
                return
            end
        
            it = 0; STOP = 1.0;
            while STOP > obj.tol0 && it < obj.itMax
                % Update iteration
                it = it + 1;
                
                % Get the residual of f, its derivative with temperature, and the
                % relative value of the residual
                [f0, fprime0, frel, molesGuess] = getRatioNewton(obj, mix1, mix2, attributeName, x0, molesGuess);
                
                % Compute pseudo-solution
                x = abs(x0 - f0 / fprime0);
                
                % Re-estimation of first derivative
                f0_2 = getRatioNewton(obj, mix1, mix2, attributeName, x, molesGuess);
                
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
            
            % Debug
            % debug_plot_error(it, aux_STOP);

            if STOP > obj.tol0
                fprintf('\n***********************************************************\n')
                fprintf('Newton method not converged\nCalling Steffensen-Aitken root finding algorithm\n')
                x0 = regulaGuess(obj, mix1, mix2, attributeName);
                [x, STOP] = obj.steff(mix1, mix2, attributeName, x0, []);
                return
            end

            printError(obj, it, x, STOP);    
        end

        function [x, STOP, molesGuess] = steff(obj, mix1, mix2, attributeName, x0, molesGuess)
            % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the Steffenson-Aitken method
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     attributeName (char): Attribute name of the problem type
            %     x0 (float): Initial guess for temperature [K]
            %     molesGuess (float): Initial guess for moles in the final mixture
            %
            % Returns:
            %     Tuple containing
            %
            %     * x (float): Temperature at equilibrium [K]
            %     * STOP (float): Relative error [-] 
            %     * molesGuess (float): Updated guess for moles in the final mixture
            
            % Check if the problem type is 'TP' or 'TV' which require no iteration
            if any(strcmpi(obj.problemType, {'TP', 'TV'}))
                x = x0;
                STOP = 0;
                return
            end
            
            % Initialization
            it = 0; STOP = 1.0;
            
            while STOP > obj.tol0 && it < obj.itMax
                % Update iteration
                it = it + 1;

                % Compute fixed point
                [g, g_rel]= getGpoint(obj, mix1, mix2, attributeName, x0, molesGuess);
                fx = abs(g - x0);

                % Compute auxiliary fixed point
                g_aux  = getGpoint(obj, mix1, mix2, attributeName, fx, molesGuess);
                fx2 = abs(g_aux - fx);

                % Compute solution
                if abs(fx2 - 2*fx + x0) > obj.tol0
                    x = obj.getPointAitken(x0, [fx, fx2]);
                else
                    x = fx;
                end
                
                % Compute stop criteria
                STOP = max(abs((x - fx) / x), abs(g_rel));

                % Update solution
                x0 = x;
            end
        
            printError(obj, it, x, STOP);
        end

        function x0 = regulaGuess(obj, mix1, mix2, attributeName)
            % Find a estimate of the temperature for the set chemical equilibrium
            % transformation using the regula falsi method
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     attributeName (char): Attribute name of the problem type
            %
            % Returns:
            %     x0 (float): Initial guess for temperature [K]
            
            % Initialization
            molesGuess = [];

            % Define branch
            x_l = obj.root_T0_l;
            x_r = obj.root_T0_r;
            
            % Compute f(x) = f2(x) - f1 at the branch limits
            g_l = obj.getGpoint(mix1, mix2, attributeName, x_l, molesGuess);
            g_r = obj.getGpoint(mix1, mix2, attributeName, x_r, molesGuess);
            
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

        function [gpoint, gpointRelative] = getGpoint(obj, mix1, mix2, attributeName, x0, molesGuess)
            % Get fixed point of a function based on the chemical transformation
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): Properties of the initial mixture
            %     mix2 (Mixture): Properties of the final mixture
            %     attributeName (str):    Fieldname in Problem Description (PD)
            %     x0 (float):     Guess temperature [K]
            %
            % Returns:
            %     Tuple containing
            %
            %     * gpoint (float): Fixed point of the function [J] (HP, EV) or [J/K] (SP, SV)
            %     * gpointRelative (float): Fixed relative point of the function [J] (HP, EV) or [J/K] (SP, SV)
            
            try
                % Compute TP problem
                obj.equilibrateT(mix1, mix2, x0, molesGuess);

                % Compute f(x) = f2(x) - f1 = 0
                gpoint = (mix2.(attributeName) - mix1.(attributeName));

                % Compute f(x) / f2(x)
                gpointRelative = gpoint / (mix2.(attributeName));
            catch
                gpoint = NaN;
                gpointRelative = NaN;
            end

        end

        function printError(obj, it, T, STOP)
            % Print error of the method if the number of iterations is greater than maximum iterations allowed
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
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