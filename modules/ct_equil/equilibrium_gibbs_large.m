function [N0, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions] = equilibrium_gibbs_large(self, pP, TP, mix1, guess_moles)
    % Obtain equilibrium composition [moles] for the given temperature [K] and pressure [bar].
    % The code stems from the minimization of the free energy of the system by using Lagrange
    % multipliers combined with a Newton-Raphson method, upon condition that initial gas
    % properties are defined by temperature and pressure.
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     pP (float): Pressure [bar]
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

    % Generalized Gibbs minimization method

    % Definitions
    N0 = self.C.N0.value;  % Composition matrix [ni, FLAG_CONDENSED_i]
    A0 = self.C.A0.value;  % Stoichiometric matrix [a_ij]
    R0TP = self.C.R0 * TP; % [J/mol]

    % Initialization
    NatomE = mix1.NatomE_react';
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
    NS0 = NS + 1;
    ind_nswt_0 = ind_nswt;
    ind_swt_0 = ind_swt;
    ind_swt = [];
    ind = [ind_nswt, ind_swt];
    NS = length(ind);

    % Initialize species vector N0
    [N0, NP] = initialize_moles(N0, NP, ind, NG, guess_moles);

    % Initialize vectors g0 zeros
    g0 = N0(:, 1);

    % Standard Gibbs free energy [J/mol]
    g0([ind_nswt_0, ind_swt_0]) = set_g0(self.S.LS([ind_nswt_0, ind_swt_0]), TP, self.DB);

    % Dimensionless chemical potential
    muRT_0 = g0/R0TP;
    muRT = muRT_0;

    % Construction of part of matrix J
    [J1, NS0] = update_matrix_J1(A0, [], NG, NS, NS0, ind, ind_elem);
    J22 = zeros(NE + 1);
    A0_T = A0';

    % Solve system
    x = equilibrium_loop;

    % Compute chemical equilibrium with condensed species
    x = equilibrium_loop_condensed(x);

    % Compute thermodynamic derivates
    [dNi_T, dN_T] = equilibrium_dT_large(self, N0, TP, A0, NG, NS, NE, ind, ind_nswt, ind_swt, ind_elem);
    [dNi_p, dN_p] = equilibrium_dp_large(self, N0, A0, NG, NS, NE, ind, ind_nswt, ind_swt, ind_elem);

    % NESTED FUNCTION
    function x = equilibrium_loop
        % Calculate composition at chemical equilibrium
        
        % Initialization
        it = 0; counter_errors = 0;
        itMax = self.TN.itMax_gibbs;
        STOP = 1.0;
        
        % Calculations
        while STOP > self.TN.tol_gibbs && it < itMax
            it = it + 1;
            % Chemical potential
            muRT(ind_nswt) =  muRT_0(ind_nswt) + log(N0(ind_nswt, 1) / NP) + log(pP);

            % Construction of matrix A (Jacobian)
            J = update_matrix_J(A0_T, J1, J22, N0, NP, ind, ind_nswt, ind_swt, ind_elem, NG);

            % Construction of vector b            
            b = update_vector_b(A0, N0, NP, NatomE, self.E.ind_E, ind_ions, ind, ind_nswt, ind_elem, muRT);

            % Solve linear system J*x = b
            x = J\b;

            % Check singular matrix
            if any(isnan(x)) || any(isinf(x))
                
                % Update temp indeces
                ind_nswt = ind_nswt_0;
                
                if FLAG_CONDENSED
                    ind_swt = ind_swt_0;
                end

                % Reset removed species to 1e-6 to try the avoid singular matrix
                N0( N0([ind_nswt, ind_swt], 1) < self.TN.tolN, 1) = 1e-6;

                if counter_errors > 2
                    x = NaN;
                    return
                end

                counter_errors = counter_errors + 1;
                continue
            end

            % Calculate correction factor
            delta = relax_factor(NP, N0(ind, 1), x(1:NS), x(end), NG);
            
            % Apply correction
            N0(ind_nswt, 1) = N0(ind_nswt, 1) .* exp(delta * x(1:NG));
            N0(ind_swt, 1) = N0(ind_swt, 1) + delta * x(NG+1:NS);
            NP = NP * exp(delta * x(end));
            
            % Compute STOP criteria
            STOP = compute_STOP(NP, x(end), N0(ind, 1), x(1:NS), NG, A0(ind, ind_elem), NatomE, max_NatomE, self.TN.tolE);
            
            % Check for negative condensed species
            if sum(N0(ind_swt, 1) < 0) 
                NS0 = NS0 + 1;
            end

            % Update temp values in order to remove species with moles < tolerance
            [ind, ind_swt, ind_nswt, ind_ions, NG, NS, N0] = update_temp(N0, ind, ind_swt, ind_nswt, ind_ions, NP, SIZE);
            
            % Update matrix J
            [J1, NS0] = update_matrix_J1(A0, J1, NG, NS, NS0, ind, ind_elem);

            % Debug            
            % aux_delta(it) = min(delta);
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
            J22 = zeros(length(ind_elem) + 1);
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
            [ind_swt_add, FLAG_CONDENSED] = check_condensed_species(A0, x(NS+1:end-1), W, ind_swt_check, muRT);
            
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
            [J1, NS0] = update_matrix_J1(A0, J1, NG, NS, NS + 1, ind, ind_elem);

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
    % Initialize composition from a previous calculation or using an
    % uniform distribution [mol]

    if isempty(guess_moles)
        N0(ind_nswt, 1) = NP/NG;
    else
        N0(ind_nswt, 1) = guess_moles(ind_nswt);
        NP = sum(guess_moles(ind_nswt));
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
    NatomE = NatomE;
    NatomE(ind_E) = 1;

    % Get flag of elements to be removed from stoichiometrix matrix
    FLAG_REMOVE_ELEMENTS = NatomE' <= tol;

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
    
    % Get indeces
    ind_elem = 1:length(NatomE);
    ind_nswt = S.ind_nswt;
    ind_swt = S.ind_swt;
    ind_ions = S.ind_nswt(S.ind_ions);
    ind = [ind_nswt, ind_swt];
    ind_cryogenic = S.ind_cryogenic;
    [ind, ind_swt] = check_cryogenic(ind, ind_swt, ind_cryogenic);
    
    % Get lengths
    NE = length(ind_elem);
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

function [J1, NS0] = update_matrix_J1(A0, J1, NG, NS, NS0, ind, ind_elem)
    % Update submatrix J1
    if NS < NS0
        J11 = eye(NS);
        J11(NG+1:end, NG+1:end) = 0;
        J12 = -[A0(ind, ind_elem), [ones(NG, 1); zeros(NS-NG, 1)]];
        J1 = [J11, J12];
        NS0 = NS;
    end

end

function [ind, ind_swt] = check_cryogenic(ind, ind_swt, ind_cryogenic)
    % Remove cryogenic species from calculations
    for i = 1:length(ind_cryogenic)
        ind(ind == ind_cryogenic(i)) = [];
        ind_swt(ind_swt == ind_cryogenic(i)) = [];
    end

end

function J2 = update_matrix_J2(A0_T, J22, N0, NP, ind, ind_nswt, ind_swt, ind_elem, NG)
    % Update submatrix J2
    J20 = N0(ind, 1)';
    J20(NG+1:end) = 0;
    J21 = [[N0(ind_nswt, 1)' .* A0_T(ind_elem, ind_nswt), A0_T(ind_elem, ind_swt)]; J20];
    J22(end, end) = -NP;
    J2 = [J21, J22];
end

function J = update_matrix_J(A0_T, J1, J22, N0, NP, ind, ind_nswt, ind_swt, ind_elem, NG)
    % Update matrix J (Jacobian)
    J2 = update_matrix_J2(A0_T, J22, N0, NP, ind, ind_nswt, ind_swt, ind_elem, NG);
    J = [J1; J2];
end

function b = update_vector_b(A0, N0, NP, NatomE, ind_E, ind_ions, ind, ind_nswt, ind_elem, muRT)
    % Update coefficient vector b
    bi = (sum(N0(ind, 1) .* A0(ind, ind_elem)))';

    if any(ind_ions)
        bi(ind_elem == ind_E) = 0;
    end

    NP_0 = NP - sum(N0(ind_nswt, 1));
    b = [-muRT(ind); NatomE - bi; NP_0];
end

function delta = relax_factor(NP, ni, eta, delta_ln_NP, NG)
    % Compute relaxation factor
    FLAG = eta(1:NG) > 0;
    FLAG_MINOR = ni(1:NG) / NP <= 1e-8 & FLAG;
    delta1 = 2./max(5*abs(delta_ln_NP), abs(eta(FLAG)));
    delta2 = min(abs((-log(ni(FLAG_MINOR)/NP) - 9.2103404) ./ (eta(FLAG_MINOR) - delta_ln_NP)));
    delta = min([1; delta1; delta2]);
end

function [STOP, deltaN1, deltaN2, deltab] = compute_STOP(NP, deltaNP, N0, deltaN0, NG, A0, NatomE, max_NatomE, tolE)
    % Compute stop criteria
    NPi = sum(N0);
    deltaN1 = N0 .* abs(deltaN0) / NPi;
    deltaN1(NG + 1:end) = abs(deltaN0(NG + 1:end)) / NPi;
    deltaN2 = NP * abs(deltaNP) / NPi;
    deltab = abs(NatomE - sum(N0 .* A0, 1)') / max_NatomE;
    deltab = max(deltab(NatomE > tolE));
    STOP = max(max(max(deltaN1), deltaN2), deltab);
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
        N0_ions = log(N0(ind_ions, 1)) + A0_ions * delta_ions;
        % Apply antilog
        N0_ions = exp(N0_ions);
        % Assign values to original vector
        N0(ind_ions, 1) = N0_ions;
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [delta_ions, ~] = ions_factor(N0, A0, ind_E, ind_nswt, ind_ions);
        % Compute STOP criteria
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