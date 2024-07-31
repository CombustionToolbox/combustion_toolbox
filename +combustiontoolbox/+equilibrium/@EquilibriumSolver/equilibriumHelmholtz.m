function [N, dNi_T, dN_T, dNi_p, dN_p, index, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, v, T, mix, molesGuess)
    % Obtain equilibrium composition [moles] for the given temperature [K] and volume [m3].
    % The code stems from the minimization of the free energy of the system by using Lagrange
    % multipliers combined with a Newton-Raphson method, upon condition that initial gas
    % properties are defined by temperature and volume.
    %
    % The algorithm implemented take advantage of the sparseness of the
    % upper left submatrix obtaining a matrix J of size NE + NS - NG. 
    %
    % This function is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311 and in Leal, A. M., Kulik, D. A.,
    % Kosakowski, G., & Saar, M. O. (2016). Computational methods for reactive transport
    % modeling: An extended law of mass-action, xLMA, method for multiphase equilibrium
    % calculations. Advances in Water Resources, 96, 405-422.
    %
    % Args:
    %     obj (EquilibriumSolver): Equilibrium solver object
    %     system (ChemicalSystem): Chemical system object
    %     v (float): Volume [m3]
    %     T (float): Temperature [K]
    %     mix1 (Mixture): Properties of the initial mixture
    %     molesGuess (float): Mixture composition [mol] of a previous computation
    %
    % Returns:
    %     Tuple containing
    %
    %     * N (float): Mixture composition [moles]
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float): Thermodynamic derivative of the moles of the mixture respect to temperature
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float): Thermodynamic derivative of the moles of the mixture respect to pressure
    %     * index (float): List of chemical species indices
    %     * STOP (float): Relative error in moles of species [-] 
    %     * STOP_ions (float): Relative error in moles of ionized species [-]
    %     * h0 (float): Molar enthalpy [J/mol]
    %
    % Examples:
    %     * N = EquilibriumSolver().equilibriumHelmholtz(system, 1, 3000, mix, [])
    %     * [N, dNi_T, dN_T, dNi_p, dN_p, index, STOP, STOP_ions] = EquilibriumSolver().equilibriumHelmholtz(system, 1, 3000, mix, [])

    % Generalized Helmholtz minimization method (reduced)
    
    % Constants
    R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]

    % Definitions
    % CHECK: ERROR IN N(:, 2) FOR CONDENSED SPECIES.
    N = system.propertyVector;       % Composition matrix [moles_i, phase_i]
    A0 = system.stoichiometricMatrix;  % Stoichiometric matrix [a_ij]
    RT = R0 * T;                       % [J/mol]
    delta0 = 0.9999;
    tau0RT = obj.tolTau;
    opts.SYM = true; % Options linsolve method: real symmetric

    % Initialization
    NatomE = mix.natomElementsReact;
    max_NatomE = max(NatomE);
    NP = 0.1;
    SIZE = -log(obj.tolMoles);
    FLAG_CONDENSED = false;
    STOP_ions = 0;

    % Set moles from molesGuess (if it was given) to 1e-6 to avoid singular matrix
    molesGuess(molesGuess < obj.tolMolesGuess) = obj.tolMolesGuess;
    
    % Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
    % for the sum of elements whose value is <= tolMoles
    [A0, indexRemoveSpecies, ind_E, NatomE] = obj.removeElements(NatomE, A0, system.ind_E, obj.tolMoles);
    
    % List of indices with nonzero values
    [index, indexGas, indexCondensed, indexIons, indexElements, NE, NG, NS] = obj.tempValues(system, NatomE);
    
    % Update temp values
    if ~isempty(indexRemoveSpecies)
        [index, indexCondensed, indexGas, indexIons, NG, NS] = obj.updateTemp(N, indexRemoveSpecies, indexCondensed, indexGas, indexIons, NP, NG, NS, SIZE);
    end

    % Remove condensed species with temperature out of bounds
    indexCondensed = check_temperature_range(system, T, indexCondensed, NS - NG, false);

    % Remove gas species with temperature out of bounds
    [indexGas, NG] = check_temperature_range(obj, T, indexGas, NG, obj.FLAG_EXTRAPOLATE);

    % First, compute chemical equilibrium with only gaseous species
    indexGas_0 = indexGas;
    indexCondensed_0 = indexCondensed;
    index0 = [indexGas_0, indexCondensed_0];
    indexCondensed = [];
    index = [indexGas, indexCondensed];
    NS = length(index);
    
    % Initialize vectors g0 (molar Gibbs energy) and h0 (molar enthalpy) with zeros
    g0 = system.propertyVector;
    h0 = system.propertyVector;
    
    % Molar Gibbs energy [J/mol]
    g0([indexGas_0, indexCondensed_0]) = set_g0(system.listSpecies([indexGas_0, indexCondensed_0]), T, system.species);
    
    % Dimensionless chemical potential
    muRT = g0/RT;
    
    % Construction of part of matrix J
    A0_T = A0';
    
    % Initialize composition matrix N [mol, FLAG_CONDENSED]
    [N, NP] = obj.equilibriumGuess(N, NP, A0_T(indexElements, index0), muRT(index0), NatomE, index0, indexGas_0, indexIons, NG, molesGuess);

    % Initialization 
    psi_j = system.propertyVector;
    tauRT = tau0RT .* min(NatomE);

    % Solve system
    x = equilibriumLoop;
    
    % Compute chemical equilibrium with condensed species
    x = equilibriumLoopCondensed(x);

    % Update matrix J (jacobian) to compute the thermodynamic derivatives
    J = update_matrix_J(A0_T(indexElements, :), N, indexGas, indexCondensed, psi_j);
    temp_zero = zeros(NS - NG + 1, 1);
    J12_2 = [sum(A0_T(indexElements, indexGas) .* N(indexGas)', 2); temp_zero(1:end-1)];
    J = [J, J12_2; J12_2', 0];

    % Molar enthalpy [J/mol]
    h0(index) = set_h0(system.listSpecies(index), T, system.species);
    
    % Dimensionless enthalpy
    H0RT = h0 / RT;

    % Compute thermodynamic derivates
    [dNi_T, dN_T, dNi_p, dN_p] = obj.equilibriumDerivatives(J, N, A0, NE, indexGas, indexCondensed, indexElements, H0RT);

    % NESTED FUNCTION
    function x = equilibriumLoop
        % Calculate composition at chemical equilibrium

        % Initialization
        it = 0; counter_errors = 0;
        itMax = obj.itMaxGibbs;
        STOP = 1.0;
        FLAG_UNSTABLE = false;
        delta_j0 = ones(NS - NG, 1);

        % Calculations
        while STOP > obj.tolGibbs && it < itMax
            it = it + 1;
            % Chemical potential
            muRT(indexGas) =  g0(indexGas) / RT + log(N(indexGas) * RT / v * 1e-5);
            
            % Compute total number of moles
            NP = sum(N(indexGas));
            
            % Construction of matrix J
            J = update_matrix_J(A0_T, N, indexGas, indexCondensed, psi_j);
            
            % Construction of vector b      
            b = update_vector_b(A0, N, NatomE, ind_E, index, indexGas, indexCondensed, indexIons, muRT, tauRT);
            
            % Solve the linear system J*x = b
            [x, ~] = linsolve(J, b, opts);
            
            % Check singular matrix
            if any(isnan(x) | isinf(x))

                % Update temp indeces
                indexGas = indexGas_0;
                indexCondensed = indexCondensed_0;
                index = [indexGas, indexCondensed];

                NG = length(indexGas);
                NS = length(index);

                % Reset removed species to 1e-6 to try the avoid singular matrix
                N( N(index) < obj.tolMoles ) = 1e-6;
                psi_j(indexCondensed) = tauRT ./ N(indexCondensed);

                if counter_errors > 2
                    x = NaN;
                    return
                end

                counter_errors = counter_errors + 1;
                continue
            end
            
            % Extract solution
            pi_i = x(1:NE);
            Delta_nj = x(NE+1:end);
            
            % Compute correction moles of gases
            Delta_ln_nj = update_Delta_ln_nj(A0, pi_i, muRT, indexGas);
            
            % Calculate correction factor
            delta = obj.relaxFactor(NP, N(index), [Delta_ln_nj; Delta_nj], 0, NG);
            
            % Apply correction
            N(indexGas) = N(indexGas) .* exp(delta * Delta_ln_nj);

            % Apply correction condensed species
            if NS - NG > 0
                delta_j = delta_j0;
                FLAG_DELTA = N(indexCondensed) + Delta_nj < 0;
                delta_j(FLAG_DELTA) = -delta0 * N(indexCondensed(FLAG_DELTA)) ./ Delta_nj(FLAG_DELTA);
                N(indexCondensed) = N(indexCondensed) + min(delta_j) .* Delta_nj;

                delta_j = delta_j0;
                Delta_psi_j = (tauRT - psi_j(indexCondensed) .* Delta_nj) ./ N(indexCondensed) - psi_j(indexCondensed);
                FLAG_DELTA = psi_j(indexCondensed) + Delta_psi_j < 0;
                delta_j(FLAG_DELTA) = -delta0 * psi_j(indexCondensed(FLAG_DELTA)) ./ Delta_psi_j(FLAG_DELTA);
                psi_j(indexCondensed) = psi_j(indexCondensed) + min(delta_j) .* Delta_psi_j;
                
                Omega_pi = exp(-psi_j(indexCondensed));
                FLAG_UNSTABLE = (N(indexCondensed) / NP < exp(-SIZE)) | (abs(log10(Omega_pi)) > 1e-2);
                N(indexCondensed(FLAG_UNSTABLE)) = 0;
            end
            
            % Compute STOP criteria
            STOP = compute_STOP(N(index), [Delta_ln_nj; Delta_nj], NG, A0(index, :), NatomE, max_NatomE, obj.tolE);
            
            % Check for negative condensed species
            FLAG_NEGATIVE = N(indexCondensed) < 0;

            % Update temp values in order to remove species with moles < tolerance
            [index, indexCondensed, indexGas, indexIons, NG, NS, N] = obj.updateTemp(N, index, indexCondensed, indexGas, indexIons, NP, NG, NS, SIZE);
            
            % Update delta_j0
            if sum(FLAG_NEGATIVE)
                FLAG_UNSTABLE(:) = false;
                delta_j0 = ones(NS - NG, 1);
            end

            % Debug 
            % aux_delta(it) = delta;
            % aux_STOP(it) = STOP;
        end

        % Check convergence of charge balance (ionized species)
        [N, STOP_ions, FLAG_ION] = equilibriumCheckIons(obj, N, A0, ind_E, indexGas, indexIons);
        
        % Additional checks in case there are ions in the mixture
        if ~FLAG_ION
            return
        end
        
        % Check that there is at least one species with n_i > tolerance 
        if any(N(indexIons) > obj.tolMoles)
            return
        end
        
        % Remove ionized species that do not satisfy n_i > tolerance
        [index, indexCondensed, indexGas, indexIons, NG, NS] = obj.updateTemp(N, index, indexCondensed, indexGas, indexIons, NP, NG, NS, SIZE);
        
        % If none of the ionized species satisfy n_i > tolerance, remove
        % electron "element" from the stoichiometric matrix
        if ~isempty(indexIons)
            return
        end
        
        % Remove element E from matrix
        indexElements(ind_E) = [];
        NE = NE - 1;

        % Debug
        % debug_plot_error(it, aux_STOP, aux_delta);
    end

    function x = equilibriumLoopCondensed(x)
        % Calculate composition at chemical equilibrium with condensed species

        if isempty(indexCondensed_0)
            return
        end

        % Update list possible gaseous species (in case singular matrix)
        indexGas_0 = indexGas;
        
        % Set list with indeces of the condensed species to be checked
        indexCondensed_check = indexCondensed_0;

        % Get molecular weight species [kg/mol]
        W = system.propertyVector;
        W(indexCondensed_check) = set_prop_DB(system.listSpecies(indexCondensed_check), 'W', system.species);

        % Definitions
        NC_max = NE - 1;
        FLAG_ALL = false;  % Include all the condensed species at once
        FLAG_ONE = false;  % Include only the condensed species that satisfies the vapour pressure test and gives the most negative value of dL_dnj
        FLAG_RULE = false; % Include only up to NC_max condensed species that satisfies the vapour pressure test with and gives the most negative values of dL_dnj

        % Initialization
        j = 0;
        while indexCondensed_check
            % Update iteration
            j = j + 1;

            % Check Gibbs phase rule
            if length(indexCondensed) > NC_max
                break;
            end

            % Check condensed species
            [indexCondensed_add, FLAG_CONDENSED, ~] = obj.equilibriumCheckCondensed(A0(:, indexElements), x(1:NE), W(indexCondensed_check), indexCondensed_check, muRT, NC_max, FLAG_ONE, FLAG_RULE);
            
            if ~FLAG_CONDENSED
                break
            end

            NC_add = length(indexCondensed_add);

            % Update indeces
            if FLAG_ALL
                indexCondensed = indexCondensed_check;
                indexCondensed_check = [];
            else
                if FLAG_ONE
                    indexCondensed_check(indexCondensed_check == indexCondensed_add) = [];
                else
                    indexCondensed_check(ismember(indexCondensed_check, indexCondensed_add)) = [];
                end
                
                indexCondensed = unique([indexCondensed, indexCondensed_add]);
            end

            index = [indexGas, indexCondensed];

            % Initialization
            STOP = 1;

            % Update lenght
            NS = length(index);

            % Save backup
            N_backup = N;

            % Check if there are non initialized condensed species
            N(indexCondensed_add( N(indexCondensed_add) == 0) ) = 1e-5;

            % Initialize Lagrange multiplier vector psi
            psi_j(indexCondensed_add) = 1e-15 ./ N(indexCondensed_add);

            % Compute chemical equilibrium considering condensed species
            x0 = equilibriumLoop;

            % Update solution vector
            if ~isnan(x0(1))
                x = x0;
                indexGas_0 = indexGas;
                continue
            end

            % Singular matrix: remove last added condensed species
            indexGas = indexGas_0;
            N = N_backup;
            N(indexCondensed(1:end-1)) = 1;
            N(indexCondensed(end)) = -1;
            [~, indexCondensed, indexGas, indexIons, NG, NS] = obj.updateTemp(N, index, indexCondensed, indexGas, indexIons, NP, NG, NS, SIZE);
            N(indexCondensed(1:end-1)) = 0;
            indexCondensed_check = indexCondensed;
        end

        % Check if there were species not considered
        [~, FLAG_CONDENSED, dL_dnj] = obj.equilibriumCheckCondensed(A0(:, indexElements), x(1:NE), W(indexCondensed_0), indexCondensed_0, muRT, NC_max, FLAG_ONE, FLAG_RULE);

        % Recompute if there are condensed species that may appear at chemical equilibrium
        if FLAG_CONDENSED && any(abs(dL_dnj) > 1e-4)
            x = equilibriumLoopCondensed(x);
        end

    end

end

% SUB-PASS FUNCTIONS
function STOP = compute_STOP(N, deltaN, NG, A0, NatomE, max_NatomE, tolE)
    % Compute stop criteria
    NPi = sum(N);
    deltaN1 = N .* abs(deltaN) / NPi;
    deltaN1(NG + 1:end) = abs(deltaN(NG + 1:end)) / NPi;
    deltab = abs(NatomE - sum(N .* A0, 1));
    deltab = max(deltab(NatomE > max_NatomE * tolE));
    STOP = max(max(deltaN1), deltab);
end

function J11 = update_matrix_J11(A0_T, N, indexGas)
    % Compute submatrix J11
    J11 = A0_T(:, indexGas) * (A0_T(:, indexGas) .* N(indexGas)')';
end

function J12 = update_matrix_J12(A0_T, indexCondensed)
    % Compute submatrix J12
    J12 = A0_T(:, indexCondensed);
end

function J = update_matrix_J(A0_T, N, indexGas, indexCondensed, psi_j)
    % Compute matrix J
    J11 = update_matrix_J11(A0_T, N, indexGas);
    J12 = update_matrix_J12(A0_T, indexCondensed);
    J22 = - diag(psi_j(indexCondensed) ./ N(indexCondensed));
    J = [J11, J12; J12', J22];
end

function b = update_vector_b(A0, N, NatomE, ind_E, index, indexGas, indexCondensed, indexIons, muRT, tauRT)
    % Compute vector b
    bi = N(index)' * A0(index, :);

    if any(indexIons)
        bi(ind_E) = NatomE(ind_E);
    end
    
    b1 = (NatomE - bi + sum(A0(indexGas, :) .* N(indexGas) .* muRT(indexGas)))';
    b2 = muRT(indexCondensed) - tauRT ./ N(indexCondensed);
    
    b = [b1; b2];
end

function Delta_ln_nj = update_Delta_ln_nj(A0, pi_i, muRT, indexGas)
    % Compute correction moles of gases
    Delta_ln_nj = sum(A0(indexGas, :)' .* pi_i, 1)' - muRT(indexGas);
end