function [N0, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions] = equilibrium_helmholtz(self, vP, TP, mix1, guess_moles)
    % Obtain equilibrium composition [moles] for the given temperature [K] and volume [m3].
    % The code stems from the minimization of the free energy of the system by using Lagrange
    % multipliers combined with a Newton-Raphson method, upon condition that initial gas
    % properties are defined by temperature and volume.
    %
    % The algorithm implemented take advantage of the sparseness of the
    % upper left submatrix obtaining a matrix J of size NE + NS - NG. 
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     vP (float): Volume [m3]
    %     TP (float): Temperature [K]
    %     mix1 (struct): Properties of the initial mixture
    %     guess_moles (float): mixture composition [mol] of a previous computation
    %
    % Returns:
    %     Tuple containing
    %
    %     * N0 (float): Equilibrium composition [moles] for the given temperature [K] and pressure [bar]
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float): Thermodynamic derivative of the moles of the mixture respect to temperature
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float): Thermodynamic derivative of the moles of the mixture respect to pressure
    %     * ind (float): List of chemical species indices
    %     * STOP (float): Relative error in moles of species [-] 
    %     * STOP_ions (float): Relative error in moles of ionized species [-] 

    % Generalized Helmholtz minimization method (reduced)
    
    % Definitions
    N0 = self.C.N0.value;  % Composition matrix [ni, FLAG_CONDENSED_i]
    A0 = self.C.A0.value;  % Stoichiometric matrix [a_ij]
    R0TP = self.C.R0 * TP; % [J/(mol)]

    % Initialization
    NatomE = mix1.NatomE;
    max_NatomE = max(NatomE);
    NP = 0.1;
    SIZE = -log(self.TN.tolN);
    FLAG_CONDENSED = false;

    % Set moles from guess_moles (if it was given) to 1e-6 to avoid singular matrix
    guess_moles(guess_moles < self.TN.tolN_guess) = self.TN.tolN_guess;

    % Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
    % for the sum of elements whose value is <= tolN
    [A0, ind_remove_species, self.E.ind_E, NatomE] = remove_elements(NatomE, A0, self.E.ind_E, self.TN.tolN);
    
    % List of indices with nonzero values
    [ind, ind_nswt, ind_swt, ind_ions, ind_elem, NE, NG, NS] = temp_values(self.S, NatomE);
    
    % Update temp values
    if ~isempty(ind_remove_species)
        [ind, ind_swt, ind_nswt, ind_ions, NG] = update_temp(N0, ind_remove_species, ind_swt, ind_nswt, ind_ions, NP, SIZE);
    end

    % Remove gas species with temperature out of bounds
    for i = length(ind_nswt):-1:1
        species = self.S.LS{ind_nswt(i)};
        if TP < self.DB.(species).T(1)
            ind_nswt(i) = [];
        end
    end

    % Remove condensed species with temperature out of bounds
    % fprintf('T = %.2f\n', TP);
    for i = length(ind_swt):-1:1
        species = self.S.LS{ind_swt(i)};
        if TP < self.DB.(species).T(1) || TP > self.DB.(species).T(end)
            ind_swt(i) = [];
        end
    end

    % First, compute chemical equilibrium with only gaseous species
    ind_nswt_0 = ind_nswt;
    ind_swt_0 = ind_swt;
    ind_swt = [];
    ind = [ind_nswt, ind_swt];
    NS = length(ind);

    % Initialize vectors g0 and h0 with zeros
    g0 = N0(:, 1);
    h0 = N0(:, 1);

    % Initialize composition matrix N0 [mol, FLAG_CONDENSED]    
    [N0, NP] = initialize_moles(N0, NP, ind_nswt, NG, guess_moles);
    
    % Standard Gibbs free energy [J/mol]
    g0([ind_nswt_0, ind_swt_0]) = set_g0(self.S.LS([ind_nswt_0, ind_swt_0]), TP, self.DB);
    
    % Dimensionless chemical potential
    muRT_0 = g0/R0TP;
    muRT = muRT_0;
    
    % Construction of part of matrix J
    J22 = [];
    A0_T = A0';
    
    % Solve system
    x = equilibrium_loop;
    
    % Compute chemical equilibrium with condensed species
    x = equilibrium_loop_condensed(x);

    % Update matrix J (jacobian) to compute the thermodynamic derivatives
    J = update_matrix_J(A0_T(ind_elem, :), J22, N0, ind_nswt, ind_swt, NE);
    temp_zero = zeros(NS - NG + 1, 1);
    J12_2 = [sum(A0_T(ind_elem, ind_nswt) .* N0(ind_nswt), 2); temp_zero(1:end-1)];
    J = [J, J12_2; J12_2', 0];
    % Standard-state enthalpy [J/mol]
    h0(ind) = set_h0(self.S.LS(ind), TP, self.DB);
    % Dimensionless standard-state enthalpy
    H0RT = h0 / R0TP;
    % Compute thermodynamic derivates
    [dNi_T, dN_T] = equilibrium_dT(J, N0, A0, NE, ind_nswt, ind_swt, ind_elem, H0RT);
    [dNi_p, dN_p] = equilibrium_dp(J, N0, A0, NE, ind_nswt, ind_swt, ind_elem);

    % NESTED FUNCTION
    function x = equilibrium_loop
        % Calculate composition at chemical equilibrium

        % Initialization
        it = 0;
        itMax = self.TN.itMax_gibbs;
        STOP = 1.0; STOP_ions = 0;

        % Calculations
        while STOP > self.TN.tol_gibbs && it < itMax
            it = it + 1;
            % Chemical potential
            muRT(ind_nswt) =  muRT_0(ind_nswt) + log(N0(ind_nswt, 1) * R0TP / vP * 1e-5);
            
            % Compute total number of moles
            NP = sum(N0(ind_nswt, 1));
            
            % Construction of matrix J
            J = update_matrix_J(A0_T, J22, N0, ind_nswt, ind_swt, NE);
            
            % Construction of vector b      
            b = update_vector_b(A0, N0, NatomE, ind, ind_nswt, ind_swt, muRT);
            
            % Solve of the linear system J*x = b
            x = J\b;
            
            % Check singular matrix
            if any(isnan(x)) || any(isinf(x))
                
                % Update temp indeces
                ind_nswt = ind_nswt_0;
                
                if FLAG_CONDENSED
                    ind_swt = ind_swt_0;
                end

                % Reset removed species to 1e-6 to try the avoid singular matrix
                N0( N0([ind_nswt, ind_swt], 1) < TN.tolN, 1) = 1e-6;

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
            Delta_ln_nj = update_Delta_ln_nj(A0, pi_i, muRT, ind_nswt);
            
            % Calculate correction factor
            delta = relax_factor(NP, N0(ind, 1), [Delta_ln_nj; Delta_nj], 0, NG);
            
            % Apply correction
            N0(ind_nswt, 1) = N0(ind_nswt, 1) .* exp(delta * Delta_ln_nj);
            N0(ind_swt, 1) = N0(ind_swt, 1) + delta * Delta_nj;
            
            % Compute STOP criteria
            STOP = compute_STOP(N0(ind, 1), [Delta_ln_nj; Delta_nj], NG, A0(ind, :), NatomE, max_NatomE, self.TN.tolE);
            
            % Check for negative condensed species
            FLAG_NEGATIVE = N0(ind_swt, 1) < 0;

            % Update temp values in order to remove species with moles < tolerance
            [ind, ind_swt, ind_nswt, ind_ions, NG, NS, N0] = update_temp(N0, ind, ind_swt, ind_nswt, ind_ions, NP, SIZE);
            
            % Update J22 matrix
            if sum(FLAG_NEGATIVE)
                J22 = zeros(NS - NG);
            end

            % Debug 
            % aux_delta(it) = delta;
            % aux_STOP(it) = STOP;
        end

        % Check convergence of charge balance (ionized species)
        [N0, STOP_ions] = check_convergence_ions(N0, A0, self.E.ind_E, ind_nswt, ind_ions, self.TN.tolN, self.TN.tol_pi_e, self.TN.itMax_ions);
        
        if ~any(N0(ind_ions) > self.TN.tolN) && ~isempty(N0(ind_ions))
            [ind, ind_swt, ind_nswt, ind_ions, NG, NS] = update_temp(N0, ind, ind_swt, ind_nswt, ind_ions, NP, SIZE);
            
            if ~isempty(ind_ions)
                return
            end
            
            % Remove element E from matrix
            ind_elem(self.E.ind_E) = [];
            NE = NE - 1;
        end

        % Debug
        % debug_plot_error(it, aux_STOP, aux_delta);
    end

    function x = equilibrium_loop_condensed(x)
        % Calculate composition at chemical equilibrium with condensed
        % species

        if isempty(ind_swt_0)
            return
        end

        % Set list with indeces of the condensed species to be checked
        ind_swt_check = ind_swt_0;

        % Get molecular weight species [kg/mol]
        W = set_prop_DB(self.S.LS(ind_swt_check), 'mm', self.DB) * 1e-3;

        % Initialization
        j = 0;
        while ind_swt_check
            j = j + 1;
            % Check condensed species
            [ind_swt_add, FLAG_CONDENSED] = check_condensed_species(A0, x(1:NE), W, ind_swt_check, muRT);
            
            if ~FLAG_CONDENSED
                break
            end

            % Update indeces
            ind_swt_check(ind_swt_check == ind_swt_add) = [];
            ind_swt = [ind_swt, ind_swt_add];
            ind = [ind_nswt, ind_swt];
            ind_swt_0 = ind_swt;

            % Initialization
            STOP = 1;

            % Update lenght
            NS = length(ind);

            % Update J matrix
            J22 = zeros(NS - NG);

            % Save backup
            N0_backup = N0;

            % Compute chemical equilibrium considering condensed species
            x0 = equilibrium_loop;

            % Update solution vector
            if ~isnan(x0(1))
                x = x0;
                ind_nswt_0 = ind_nswt;
                continue
            end

            % Singular matrix: remove last added condensed species
            ind_nswt = ind_nswt_0;
            N0 = N0_backup;
            [~, ind_swt, ind_nswt, ind_ions, NG, NS] = update_temp(N0, ind_swt(end), ind_swt, ind_nswt, ind_ions, NP, SIZE);
        end

    end

end

% SUB-PASS FUNCTIONS
function [N0, NP] = initialize_moles(N0, NP, ind_nswt, NG, guess_moles)
    % Initialize composition [mol]
    if isempty(guess_moles)
        N0(ind_nswt, 1) = NP / NG;
    else
        N0(ind_nswt, 1) = guess_moles(ind_nswt);
        NP = sum(guess_moles);
    end
end

function [ind_swt, FLAG_CONDENSED] = check_condensed_species(A0, pi_i, W, ind_swt, muRT)
    % Check condensed species
    
    % Initialization
    FLAG_CONDENSED = false;
    % Get length condensed species
    NC = length(ind_swt);

    for i = NC:-1:1
        % Only check if there were atoms of the species in the initial
        % mixture
        if ~sum(A0(ind_swt(i), :))
            continue
        end

        % Calculate dLdnj of the condensed species
        dL_dn(i) = (muRT(ind_swt(i)) - dot(pi_i, A0(ind_swt(i), :))) / W(i);
    end
    
    FLAG = dL_dn < 0;
    % Check if any condensed species have to be considered
    if ~sum(FLAG)
        ind_swt = [];
        return
    end
    
    % Get index of the condensed species to be added to the system
    [~, temp] = min(dL_dn);
    ind_swt = ind_swt(temp);
    % Update flag
    FLAG_CONDENSED = true;
end

function ind_remove_species = find_remove_species(A0, FLAG_REMOVE_ELEMENTS)
    % Get flag of species to be removed from stoichiometrix matrix A0
    ind_remove_species = find(sum(A0(:, FLAG_REMOVE_ELEMENTS) > 0, 2) > 0);
end

function [A0, ind_remove_species, ind_E, NatomE] = remove_elements(NatomE, A0, ind_E, tol)
    % Find zero sum elements

    % Define temporal fictitious value if there are ionized species
    temp_NatomE = NatomE;
    temp_NatomE(ind_E) = 1;

    % Get flag of elements to be removed from stoichiometrix matrix
    FLAG_REMOVE_ELEMENTS = temp_NatomE' <= tol;

    % Get the species to be removed from stoichiometrix matrix
    ind_remove_species = find_remove_species(A0, FLAG_REMOVE_ELEMENTS);

    % Update stoichiometrix matrix
    A0(:, FLAG_REMOVE_ELEMENTS) = [];

    % Set number of atoms
    NatomE(FLAG_REMOVE_ELEMENTS) = [];

    % Check position "element" electron
    if ind_E
        ind_E = ind_E - sum(FLAG_REMOVE_ELEMENTS(1:ind_E-1));
    end
end

function [ind, ind_nswt, ind_swt, ind_ions, ind_elem, NE, NG, NS] = temp_values(S, NatomE)
    % List of indices with nonzero values and lengths
    ind_elem = 1:length(NatomE);
    ind_nswt = S.ind_nswt;
    ind_swt = S.ind_swt;
    ind_ions = S.ind_nswt(S.ind_ions);
    ind_cryogenic = S.ind_cryogenic;
    ind = [ind_nswt, ind_swt];
    [ind, ind_swt] = check_cryogenic(ind, ind_swt, ind_cryogenic);

    % Update lengths
    NE = length(NatomE);
    NG = length(ind_nswt);
    NS = length(ind);
end

function [ind_swt, ind_nswt, ind_ions, n] = remove_item(n, ind, ind_swt, ind_nswt, ind_ions, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if n(i) / NP < exp(-SIZE)
            n(i) = 0;
            ind_swt(ind_swt==ind(i)) = [];
            ind_nswt(ind_nswt==ind(i)) = [];
            ind_ions(ind_ions==ind(i)) = [];
        end
        
    end

end

function [ind, ind_swt, ind_nswt, ind_ions, NG, NS, N0] = update_temp(N0, ind, ind_swt, ind_nswt, ind_ions, NP, SIZE)
    % Update temp items
    [ind_swt, ind_nswt, ind_ions, n] = remove_item(N0(ind, 1), ind, ind_swt, ind_nswt, ind_ions, NP, SIZE);
    N0(ind, 1) = n;
    ind = [ind_nswt, ind_swt];
    NG = length(ind_nswt);
    NS = length(ind);
end

function [ind, ind_swt] = check_cryogenic(ind, ind_swt, ind_cryogenic)
    % Remove cryogenic species from calculations
    for i = 1:length(ind_cryogenic)
        ind(ind == ind_cryogenic(i)) = [];
        ind_swt(ind_swt == ind_cryogenic(i)) = [];
    end

end

function delta = relax_factor(NP, ni, eta, Delta_ln_NP, NG)
    % Compute relaxation factor
    FLAG = eta(1:NG) > 0;
    FLAG_MINOR = ni(1:NG) / NP <= 1e-8 & FLAG;
    delta1 = 2./max(5*abs(Delta_ln_NP), abs(eta(FLAG)));
    delta2 = min(abs((-log(ni(FLAG_MINOR)/NP) - 9.2103404) ./ (eta(FLAG_MINOR) - Delta_ln_NP)));
    delta = min([1; delta1; delta2]);
end

function STOP = compute_STOP(N0, deltaN0, NG, A0, NatomE, max_NatomE, tolE)
    % Compute stop criteria
    NPi = sum(N0);
    deltaN1 = N0 .* abs(deltaN0) / NPi;
    deltaN1(NG + 1:end) = abs(deltaN0(NG + 1:end)) / NPi;
    deltab = abs(NatomE - sum(N0 .* A0, 1));
    deltab = max(deltab(NatomE > max_NatomE * tolE));
    STOP = max(max(deltaN1), deltab);
end

function J11 = update_matrix_J11(A0_T, N0, ind_nswt, NE)
    % Compute submatrix J11
    for k = NE:-1:1
        J11(:, k) = sum(A0_T(k, ind_nswt) .* A0_T(:, ind_nswt) .* N0(ind_nswt), 2);
    end
end

function J12 = update_matrix_J12(A0_T, ind_swt)
    % Compute submatrix J12
    J12 = A0_T(:, ind_swt);
end

function J = update_matrix_J(A0_T, J22, N0, ind_nswt, ind_swt, NE)
    % Compute matrix J
    J11 = update_matrix_J11(A0_T, N0, ind_nswt, NE);
    J12 = update_matrix_J12(A0_T, ind_swt);
    J = [J11, J12; J12', J22];
end

function b = update_vector_b(A0, N0, NatomE, ind, ind_nswt, ind_swt, muRT) 
    % Compute vector b
    bi_0 = (NatomE - N0(ind, 1)' * A0(ind, :) + sum(A0(ind_nswt, :) .* N0(ind_nswt, 1) .* muRT(ind_nswt)))';    
    b = [bi_0; muRT(ind_swt)];
end

function Delta_ln_nj = update_Delta_ln_nj(A0, pi_i, muRT, ind_nswt)
    % Compute correction moles of gases
    Delta_ln_nj = sum(A0(ind_nswt, :)' .* pi_i)' - muRT(ind_nswt);
end

function [N0, STOP] = check_convergence_ions(N0, A0, ind_E, ind_nswt, ind_ions, TOL, TOL_pi, itMax)
    % Check convergence of ionized species

    % Initialization
    STOP = 0;

    % Check if there are ionized species
    if ~any(ind_ions)
        return
    end
    
    % Get error in the electro-neutrality of the mixture
    [delta_ions, ~] = ions_factor(N0, A0, ind_E, ind_nswt, ind_ions);
    
    % Reestimate composition of ionized species
    if abs(delta_ions) > TOL_pi
        [N0, STOP] = recompute_ions(N0, A0, ind_E, ind_nswt, ind_ions, delta_ions, TOL, TOL_pi, itMax);
    end
    
end

function [N0, STOP] = recompute_ions(N0, A0, ind_E, ind_nswt, ind_ions, delta_ions, TOL, TOL_pi, itMax)
    % Reestimate composition of ionized species
    
    % Initialization
    A0_ions = A0(ind_ions, ind_E);
    STOP = 1;
    it = 0;
    % Reestimate composition of ionized species
    while STOP > TOL_pi && it < itMax
        it = it + 1;
        % Apply correction
        N0(ind_ions, 1) = N0(ind_ions, 1) .* exp(A0_ions * delta_ions);
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [delta_ions, ~] = ions_factor(N0, A0, ind_E, ind_nswt, ind_ions);
        STOP = abs(delta_ions);
    end   
    
    Xi_ions = N0(ind_ions, 1) / sum(N0(:, 1));

    % Set error to zero if molar fraction of ionized species are below
    % tolerance
    if ~any(Xi_ions > TOL)
        STOP = 0;
    end

end

function [delta, deltaN3] = ions_factor(N0, A0, ind_E, ind_nswt, ind_ions)
    % Compute relaxation factor for ionized species
    
    if ~any(ind_ions)
        delta = [];
        deltaN3 = 0;
        return
    end

    delta = -sum(A0(ind_nswt, ind_E) .* N0(ind_nswt, 1))/ ...
             sum(A0(ind_nswt, ind_E).^2 .* N0(ind_nswt, 1));
    deltaN3 = abs(sum(N0(ind_nswt, 1) .* A0(ind_nswt, ind_E)));
end