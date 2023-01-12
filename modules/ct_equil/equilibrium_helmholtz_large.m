function [N0, dNi_T, dN_T, dNi_p, dN_p, STOP, STOP_ions] = equilibrium_helmholtz_large(self, vP, TP, mix1, guess_moles)
    % Obtain equilibrium composition [moles] for the given temperature [K] and volume [m3].
    % The code stems from the minimization of the free energy of the system by using Lagrange
    % multipliers combined with a Newton-Raphson method, upon condition that initial gas
    % properties are defined by temperature and volume.
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
    %     * STOP (float): Relative error in moles of species [-] 
    %     * STOP_ions (float): Relative error in moles of ionized species [-] 

    % Generalized Helmholtz minimization method
    
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
    [temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_elem, temp_NE, temp_NG, temp_NS] = temp_values(self.S, NatomE);
    
    % Update temp values
    if ~isempty(ind_remove_species)
        [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NG] = update_temp(N0(ind_remove_species, 1), ind_remove_species, temp_ind_swt, temp_ind_nswt, temp_ind_ions, NP, SIZE);
    end
    
    % Update temp values
    temp_NS0 = temp_NS + 1;
    
    temp_ind_nswt_0 = temp_ind_nswt;
    temp_ind_swt_0 = temp_ind_swt;
    temp_ind_swt = [];
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NS = length(temp_ind);

    % Initialize species vector N0    
    [N0, NP] = initialize_moles(N0, NP, temp_ind, temp_NG, guess_moles);
    
    % Standard Gibbs free energy [J/mol]
    g0 = set_g0(self.S.LS, TP, self.DB);

    % Dimensionless chemical potential
    muRT_0 = g0/R0TP;
    muRT = muRT_0;

    % Get Molar mass [kg/mol]
    % for i = self.S.NS:-1:1
    %     W(i) = self.DB.(self.S.LS{i}).mm * 1e-3;
    % end

    % Construction of part of matrix J
    [J1, temp_NS0] = update_matrix_J1(A0, [], temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_elem);
    J22 = zeros(temp_NE);
    A0_T = A0';
    
    % Solve system
    x = equilibrium_loop;

    % Check condensed species
    [temp_ind, temp_ind_swt, FLAG] = check_condensed_species(A0, x, temp_ind, temp_ind_nswt, temp_ind_swt_0, temp_ind_elem, temp_NE, muRT);
    if FLAG
        STOP = 1;
        temp_NS = length(temp_ind);
        temp_NS0 = temp_NS + 1;
        [J1, temp_NS0] = update_matrix_J1(A0, J1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_elem);

        self.TN.itMax_gibbs = self.TN.itMax_gibbs / 2;
        equilibrium_loop;
    end

    % Compute thermodynamic derivates
    [dNi_T, dN_T] = equilibrium_dT_large(self, N0, TP, A0, temp_NG, temp_NS, temp_NE, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_elem);
    [dNi_p, dN_p] = equilibrium_dp_large(self, N0, A0, temp_NG, temp_NS, temp_NE, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_elem);

    % NESTED FUNCTION
    function x = equilibrium_loop
        % Calculate composition at chemical equilibrium
        
        % Initialization
        it = 0; counter_errors = 0;
        itMax = self.TN.itMax_gibbs;
        STOP = 1;
        
        % Calculations
        while STOP > self.TN.tol_gibbs && it < itMax
            it = it + 1;
            % Chemical potential
            muRT(temp_ind_nswt) =  muRT_0(temp_ind_nswt) + log(N0(temp_ind_nswt, 1) * R0TP / vP * 1e-5);
            
            % Compute total number of moles
            NP = sum(N0(temp_ind_nswt, 1));
            
            % Helmholtz free energy [cal/g] (debug)
            % Helmholtz(it) = dot(N0(:, 1), muRT * R0TP) / dot(N0(:, 1), W) / 4186.8;
            % fprintf('Helmholtz: %f\n', Helmholtz(it));
            
            % Construction of matrix J
            J = update_matrix_J(A0_T, J1, J22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_elem);

            % Construction of vector b            
            b = update_vector_b(A0, N0, NatomE, self.E.ind_E, temp_ind_ions, temp_ind, temp_ind_elem, muRT);
            
            % Solve linear system J*x = b
            x = J\b;

            % Check singular matrix
            if any(isnan(x)) || any(isinf(x))
                
                % Update temp indeces
                temp_ind_nswt = temp_ind_nswt_0;
                
                if FLAG_CONDENSED
                    temp_ind_swt = temp_ind_swt_0;
                end

                % Reset removed species to 1e-6 to try the avoid singular matrix
                N0( N0([temp_ind_nswt, temp_ind_swt], 1) < self.TN.tolN, 1) = 1e-6;

                if counter_errors > 2
                    x = NaN;
                    return
                end

                counter_errors = counter_errors + 1;
                continue
            end

            % Calculate correction factor
            delta = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), 0, temp_NG);
            
            % Apply correction
            N0(temp_ind_nswt, 1) = N0(temp_ind_nswt, 1) .* exp(delta * x(1:temp_NG));
            N0(temp_ind_swt, 1) = N0(temp_ind_swt, 1) + delta * x(temp_NG+1:temp_NS);
            
            % Compute STOP criteria
            STOP = compute_STOP(N0(temp_ind, 1), x(1:temp_NS), temp_NG, A0(temp_ind, temp_ind_elem), NatomE, max_NatomE, self.TN.tolE);
            
            % Update temp values in order to remove species with moles < tolerance
            [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NG, temp_NS] = update_temp(N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, NP, SIZE);
            
            % Update matrix J
            [J1, temp_NS0] = update_matrix_J1(A0, J1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_elem);
            % Debug            
            % aux_lambda(it) = min(lambda);
            % aux_STOP(it) = STOP;
        end

        % Check convergence of charge balance (ionized species)
        [N0, STOP_ions] = check_convergence_ions(N0, A0, self.E.ind_E, temp_ind_nswt, temp_ind_ions, self.TN.tolN, self.TN.tol_pi_e, self.TN.itMax_ions);
        if ~any(N0(temp_ind_ions) > self.TN.tolN) && ~isempty(N0(temp_ind_ions))
            temp_ind_elem(temp_ind_elem == self.E.ind_E) = [];
            temp_NE = temp_NE - 1;
        end

        % Debug  
        % debug_plot_error(it, aux_STOP, aux_lambda);
    end

end

% SUB-PASS FUNCTIONS
function [N0, NP] = initialize_moles(N0, NP, temp_ind_nswt, temp_NG, guess_moles)
    % Initialize composition from a previous calculation or using an
    % uniform distribution [mol]

    if isempty(guess_moles)
        N0(temp_ind_nswt, 1) = NP/temp_NG;
    else
        N0(temp_ind_nswt, 1) = guess_moles(temp_ind_nswt);
        NP = sum(guess_moles(temp_ind_nswt));
    end

end

function [temp_ind, temp_ind_swt, FLAG_CONDENSED] = check_condensed_species(A0, x, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_elem, temp_NE, muRT)
    % Check condensed species
    
    % Initialization
    FLAG_CONDENSED = false;
    
    % Check if there are condensed species
    if isempty(temp_ind_swt)
        return
    end
    
    % Get length condensed species
    NC = length(temp_ind_swt);
    
    % Initialize false vector
    temp = false(NC, 1);
    
    for i = NC:-1:1
        % Only check if there were atoms of the species in the initial
        % mixture
        if ~sum(A0(temp_ind_swt(i), temp_ind_elem))
            continue
        end

        dG_dn = muRT(temp_ind_swt(i)) - dot(x(end-temp_NE+1:end), A0(temp_ind_swt(i), temp_ind_elem));
        
        if dG_dn < 0
            temp(i) = true;
        end

    end
    
    % Check if any condensed species have to be considered
    if ~sum(temp)
        temp_ind_swt = [];
        return
    end
    
    % Update flag
    FLAG_CONDENSED = true;

    % Update indeces
    temp_ind_swt = temp_ind_swt(temp);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
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

function [temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_elem, temp_NE, temp_NG, temp_NS] = temp_values(S, NatomE)
    % List of indices with nonzero values and lengths
    
    % Get indeces
    temp_ind_elem = 1:length(NatomE);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_ind_ions = S.ind_nswt(S.ind_ions);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_ind_cryogenic = S.ind_cryogenic;
    [temp_ind, temp_ind_swt] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic);
    
    % Get lengths
    temp_NE = length(temp_ind_elem);
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end

function [temp_ind_swt, temp_ind_nswt, temp_ind_ions] = remove_item(n, ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if log(n(i)/NP) < -SIZE
            temp_ind_swt(temp_ind_swt==ind(i)) = [];
            temp_ind_nswt(temp_ind_nswt==ind(i)) = [];
            temp_ind_ions(temp_ind_ions==ind(i)) = [];   
        end

    end

end

function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NG, temp_NS] = update_temp(n, ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, NP, SIZE)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt, temp_ind_ions] = remove_item(n, ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, NP, SIZE);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end

function [J1, temp_NS0] = update_matrix_J1(A0, J1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_elem)
    % Update stoichiometric submatrix J1
    if temp_NS < temp_NS0
        J11 = eye(temp_NS);
        J11(temp_NG+1:end, temp_NG+1:end) = 0;
        J12 = -A0(temp_ind, temp_ind_elem);
        J1 = [J11, J12];
        temp_NS0 = temp_NS;
    end
end

function [temp_ind, temp_ind_swt] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic)
    % Remove cryogenic species from calculations
    for i = 1:length(temp_ind_cryogenic)
        temp_ind(temp_ind == temp_ind_cryogenic(i)) = [];
        temp_ind_swt(temp_ind_swt == temp_ind_cryogenic(i)) = [];
    end

end

function J2 = update_matrix_J2(A0_T, J22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_elem)
    % Update stoichiometric submatrix J2
    J21 = [N0(temp_ind_nswt, 1)' .* A0_T(temp_ind_elem, temp_ind_nswt), A0_T(temp_ind_elem, temp_ind_swt)];
    J2 = [J21, J22];
end

function J = update_matrix_J(A0_T, J1, J22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_elem)
    % Update stoichiometric matrix J
    J2 = update_matrix_J2(A0_T, J22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_elem);
    J = [J1; J2];
end

function b = update_vector_b(A0, N0, NatomE, ind_E, temp_ind_ions, temp_ind, temp_ind_elem, muRT)
    % Update coefficient vector b
    bi = (sum(N0(temp_ind, 1) .* A0(temp_ind, temp_ind_elem)))';
    
    if any(temp_ind_ions)
        bi(temp_ind_elem == ind_E) = 0;
    end

    b = [-muRT(temp_ind); NatomE - bi];
end

function lambda = relax_factor(NP, ni, eta, Delta_ln_NP, temp_NG)
    % Compute relaxation factor
    FLAG = eta(1:temp_NG) > 0;
    FLAG_MINOR = ni(1:temp_NG) / NP <= 1e-8 & FLAG;
    lambda1 = 2./max(5*abs(Delta_ln_NP), abs(eta(FLAG)));
    lambda2 = min(abs((-log(ni(FLAG_MINOR)/NP) - 9.2103404) ./ (eta(FLAG_MINOR) - Delta_ln_NP)));
    lambda = min([1; lambda1; lambda2]);
end

function [STOP, deltaN1, deltab] = compute_STOP(N0, deltaN0, temp_NG, A0, NatomE, max_NatomE, tolE)
    % Compute stop criteria
    NPi = sum(N0);
    deltaN1 = N0 .* abs(deltaN0) / NPi;
    deltaN1(temp_NG+1:end) = abs(deltaN0(temp_NG+1:end)) / NPi;
    deltab = abs(NatomE - sum(N0 .* A0, 1)');
    deltab = max(deltab(NatomE > max_NatomE * tolE));
    STOP = max(max(deltaN1), deltab);
end

function [N0, STOP] = check_convergence_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, TOL, TOL_pi, itMax)
    % Check convergence of ionized species
    
    % Initialization
    STOP = 0;

    % Check if there are ionized species
    if ~any(temp_ind_ions)
        return
    end
    
    % Get error in the electro-neutrality of the mixture
    [delta_ions, ~] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions);
    
    % Reestimate composition of ionized species
    if abs(delta_ions) > TOL_pi
        [N0, STOP] = recompute_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, delta_ions, TOL, TOL_pi, itMax);
    end

end

function [N0, STOP] = recompute_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, delta_ions, TOL, TOL_pi, itMax)
    % Reestimate composition of ionized species
    
    % Initialization
    A0_ions = A0(temp_ind_ions, ind_E);
    STOP = 1;
    it = 0;
    % Reestimate composition of ionized species
    while STOP > TOL_pi && it < itMax
        it = it + 1;
        % Apply correction
        N0_ions = log(N0(temp_ind_ions, 1)) + A0_ions * delta_ions;
        % Apply antilog
        N0_ions = exp(N0_ions);
        % Assign values to original vector
        N0(temp_ind_ions, 1) = N0_ions;
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [delta_ions, ~] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions);
        % Compute STOP criteria
        STOP = abs(delta_ions);
    end   
    
    Xi_ions = N0(temp_ind_ions, 1) / sum(N0(:, 1));
    
    % Set error to zero if molar fraction of ionized species are below
    % tolerance
    if ~any(Xi_ions > TOL)
        STOP = 0;
    end

end

function [delta, deltaN3] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions)
    % Compute relaxation factor for ionized species
    
    if ~any(temp_ind_ions)
        delta = [];
        deltaN3 = 0;
        return
    end

    delta = -sum(A0(temp_ind_nswt, ind_E) .* N0(temp_ind_nswt, 1))/ ...
             sum(A0(temp_ind_nswt, ind_E).^2 .* N0(temp_ind_nswt, 1));
    deltaN3 = abs(sum(N0(temp_ind_nswt, 1) .* A0(temp_ind_nswt, ind_E)));
end