function [N0, dNi_T, dN_T, dNi_p, dN_p, STOP, STOP_ions] = equilibrium_helmholtz_eos(self, vP, TP, mix1, guess_moles)
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
    %     Tuple containing:
    %
    %     * N0 (float): Equilibrium composition [moles] for the given temperature [K] and pressure [bar]
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float): Thermodynamic derivative of the moles of the mixture respect to temperature
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float): Thermodynamic derivative of the moles of the mixture respect to pressure
    %     * STOP (float): Relative error in moles of species [-] 
    %     * STOP_ions (float): Relative error in moles of ionized species [-] 

    % Generalized Helmholtz minimization method
    
    % Abbreviations ---------------------
    E = self.E;
    S = self.S;
    C = self.C;
    TN = self.TN;
    % -----------------------------------
    % Definitions
    N0 = C.N0.value;
    A0 = C.A0.value;
    R0TP = C.R0 * TP; % [J/mol]
    % Initialization
    NatomE = mix1.NatomE_react';
    max_NatomE = max(NatomE);
    NP = 0.1;
    SIZE = -log(TN.tolN);
    flag_ions_first = true;

    % Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
    % for the sum of elements whose value is <= tolN
    [A0, ind_A0_E0, E.ind_E, NatomE] = remove_elements(NatomE, A0, S.LS, E.ind_E, guess_moles, TN.tolN);
    % List of indices with nonzero values
    [temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_E, temp_NE, temp_NG, temp_NS] = temp_values(E.ind_E, S, NatomE, TN.tolN);
    % Update temp values
    [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NG] = update_temp(N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, temp_ind_E, E.ind_E, temp_ind_ions, [], NP, SIZE, flag_ions_first);

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
    g0 = set_g0(S.LS, TP, self.DB);
    % Dimensionless Chemical potential
    muRT_0 = g0/R0TP;
    muRT = muRT_0;
%     % Get Molar mass [kg/mol]
%     for i = S.NS:-1:1
%         W(i) = self.DB.(S.LS{i}).mm * 1e-3;
%     end
    % Construction of part of matrix A
    [A1, temp_NS0] = update_matrix_A1(A0, [], temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
    A22 = zeros(temp_NE);
    A0_T = A0';
    % Solve system
    x = equilibrium_loop;
    % Check condensed species
    [temp_ind, temp_ind_swt, FLAG] = check_condensed_species(A0, x, temp_ind_nswt, temp_ind_swt_0, temp_ind_E, temp_NE, muRT);
    if FLAG
        if any(isnan(x))
            % Update temp values
            [temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_E, temp_NE, temp_NG, temp_NS] = temp_values(E.ind_E, S, NatomE, TN.tolN);
            % Initialize species vector N0 
            N0(:, 1) = NP/temp_NS;
            % Construction of part of matrix A
            A22 = zeros(temp_NE + 1);
            A0_T = A0';
        end
        STOP = 1;
        temp_NS = length(temp_ind);
        temp_NS0 = temp_NS + 1;
        [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);

        TN.itMax_gibbs = TN.itMax_gibbs / 2;
        equilibrium_loop;
    end

    % Compute thermodynamic derivates
    [dNi_T, dN_T] = equilibrium_dT(self, N0, TP, A0, temp_NG, temp_NS, temp_NE, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E);
    [dNi_p, dN_p] = equilibrium_dp(self, N0, A0, temp_NG, temp_NS, temp_NE, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E);

    % NESTED FUNCTION
    function x = equilibrium_loop
        it = 0;
%         itMax = 50 + round(S.NS/2);
        itMax = TN.itMax_gibbs;
        STOP = 1;
        while STOP > TN.tolN && it < itMax
            it = it + 1;
            % Non-ideal contribution of chemical potential
            muRT_ex = mu_ex_ideal(self, N0(temp_ind_nswt, 1), TP, vP) / R0TP;
            % Chemical potential
            muRT(temp_ind_nswt) =  muRT_0(temp_ind_nswt) + log(N0(temp_ind_nswt, 1) * R0TP / vP * 1e-5) + muRT_ex;
            % Compute total number of moles
            NP = sum(N0(temp_ind_nswt, 1));
            % Helmholtz free energy [cal/g] (debug)
%             Helmholtz(it) = dot(N0(:, 1), muRT * R0TP) / dot(N0(:, 1), W) / 4186.8;
%             fprintf('Helmholtz: %f\n', Helmholtz(it));
            % Construction of matrix A
            A = update_matrix_A(A0_T, A1, A22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E);
            % Check singularity
%             A = check_singularity(A, it);
            % Construction of vector b            
            b = update_vector_b(A0, N0, NatomE, E.ind_E, temp_ind_ions, temp_ind, temp_ind_E, muRT);
            % Solve of the linear system A*x = b
            x = A\b;
            % Omit nan or inf values
            x(isnan(x) | isinf(x)) = 0;
            % Calculate correction factor
            lambda = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), 0, temp_NG);
            % Apply correction
            N0(temp_ind_nswt, 1) = exp(log(N0(temp_ind_nswt, 1)) + lambda * x(1:temp_NG));
            N0(temp_ind_swt, 1) = N0(temp_ind_swt, 1) + lambda * x(temp_NG+1:temp_NS);
            % Compute STOP criteria
            STOP = compute_STOP(N0(temp_ind, 1), x(1:temp_NS), temp_NG, A0(temp_ind, temp_ind_E), NatomE, max_NatomE, TN.tolE);
            % Update temp values in order to remove species with moles < tolerance
            [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NG, temp_NS, temp_ind_E, A22, flag_ions_first, N0] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, E.ind_E, temp_ind_ions, A22, NP, SIZE, flag_ions_first);
            % Update matrix A
            [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
            % Debug            
%             aux_lambda(it) = min(lambda);
%             aux_STOP(it) = STOP;
        end
        % Check convergence of charge balance (ionized species)
        [N0, STOP_ions] = check_convergence_ions(N0, A0, E.ind_E, temp_ind_nswt, temp_ind_ions, TN.tolN, TN.tol_pi_e, TN.itMax_ions);
        if ~any(N0(temp_ind_ions) > TN.tolN) && ~isempty(N0(temp_ind_ions))
            temp_ind_E(temp_ind_E == E.ind_E) = [];
            temp_NE = temp_NE - 1;
        end
        % Debug  
%         debug_plot_error(it, aux_STOP, aux_lambda);
    end

    function A = check_singularity(A, it)
        % Check if matrix is matrix is singular
        if it > 10*temp_NS
            if cond(A) > exp(0.5 * log(TN.tolN))
                FLAG_OLD = ~ismember(temp_ind_nswt_0, temp_ind_nswt);
                if any(FLAG_OLD)
                    temp_ind_nswt = temp_ind_nswt_0;
                    temp_ind = [temp_ind_nswt, temp_ind_swt];
                    temp_NG = length(temp_ind_nswt);
                    temp_NS = length(temp_ind);
                    N0(temp_ind_nswt(FLAG_OLD)) = TN.tolN * 1e-1;
                    [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
                    A = update_matrix_A(A0_T, A1, A22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E);
                    N0(temp_ind_nswt(FLAG_OLD)) = 0;
                end
            end
        end
    end
end

% SUB-PASS FUNCTIONS
function [N0, NP] = initialize_moles(N0, NP, temp_ind_nswt, temp_NG, guess_moles)
    if isempty(guess_moles)
        N0(temp_ind_nswt, 1) = NP/temp_NG;
    else
        N0(temp_ind_nswt, 1) = guess_moles(temp_ind_nswt);
        NP = sum(guess_moles);
    end
end

function [temp_ind, temp_ind_swt, FLAG] = check_condensed_species(A0, x, temp_ind_nswt, temp_ind_swt_0, temp_ind_E, temp_NE, muRT)
    % Check condensed species
    aux = [];
    FLAG = false;
    if any(isnan(x))
        aux = true;
    else
        for i=length(temp_ind_swt_0):-1:1
            dG_dn = muRT(temp_ind_swt_0(i)) - dot(x(end-temp_NE+1:end), A0(temp_ind_swt_0(i), temp_ind_E));
            if dG_dn < 0
                aux = [aux, i];
            end
        end
    end
    if any(aux)
        FLAG = true;
    end
    if ~isempty(temp_ind_swt_0)
        temp_ind_swt = temp_ind_swt_0(aux);
    else
        temp_ind_swt = [];
    end
    temp_ind = [temp_ind_nswt, temp_ind_swt];
end

function ind_A = find_ind_Matrix(A, bool)
    ls = find(bool>0);
    ind_A = [];
    i = 1;
    for ind = ls
        ind_A = [ind_A, find(A(:, ind) > 0)'];
        i = i + 1;
    end
end

function [A0, ind_A0_E0, ind_E, NatomE] = remove_elements(NatomE, A0, species, ind_E, guess_moles, tol)
    % Find zero sum elements

    % Ions
    temp_ind_ions = contains(species, 'minus') | contains(species, 'plus'); %S.ind_ions;
    if ~isempty(guess_moles)
        pass_ions = guess_moles(temp_ind_ions) > TN.tolN;
        temp_ind_ions = temp_ind_ions(pass_ions);
    end
    aux = NatomE;
    if any(temp_ind_ions)
        NatomE(ind_E) = 1; % temporal fictitious value
    end
    % No-ions
    bool_E0 = NatomE' <= tol;
    ind_A0_E0 = find_ind_Matrix(A0, bool_E0);
    A0 = A0(:, NatomE > tol);
    % Set number of atoms
    if any(temp_ind_ions)
        bool_range1 = aux(1:ind_E-1) > tol;
        bool_range2 = aux(ind_E+1:end) > tol;
        NatomE = aux([bool_range1; true; bool_range2]);
        ind_E = sum(bool_range1) + 1;
    else
        NatomE = aux(aux > tol);
    end
end

function [temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_E, temp_NE, temp_NG, temp_NS] = temp_values(ind_E, S, NatomE, tol)
    % List of indices with nonzero values and lengths
    FLAG_IONS = S.ind_ions;
    if any(FLAG_IONS)
        NatomE(ind_E) = 1; % temporal fictitious value
    end
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_ind_ions = S.ind_nswt(FLAG_IONS);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_ind_cryogenic = S.ind_cryogenic;
    [temp_ind, temp_ind_swt] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic);

    temp_NE = length(temp_ind_E);
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end

function [temp_ind_swt, temp_ind_nswt, temp_ind_E, temp_ind_ions, A22, flag_ions_first, N0] = remove_item(N0, n, ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, ind_E, temp_ind_ions, A22, NP, SIZE, flag_ions_first)
    % Remove species from the computed indeces list of gaseous and condupdate_tempensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if log(n(i)/NP) < -SIZE
            if N0(ind(i), 2)
                temp_ind_swt(temp_ind_swt==ind(i)) = [];
            else
                temp_ind_nswt(temp_ind_nswt==ind(i)) = [];
                temp_ind_ions(temp_ind_ions==ind(i)) = [];
                if isempty(temp_ind_ions) && flag_ions_first
                    % remove element E from matrix
                    temp_ind_E(ind_E) = [];
                    A22 = zeros(length(temp_ind_E));
                    flag_ions_first = false;
                end
            end
%             N0(ind(i), 1) = 0;
        end
    end
end

function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NG, temp_NS, temp_ind_E, A22, flag_ions_first, N0] = update_temp(N0, zip1, zip2, temp_ind_swt, temp_ind_nswt, temp_ind_E, ind_E, temp_ind_ions, A22, NP, SIZE, flag_ions_first)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt, temp_ind_E, temp_ind_ions, A22, flag_ions_first, N0] = remove_item(N0, zip1, zip2, temp_ind_swt, temp_ind_nswt, temp_ind_E, ind_E, temp_ind_ions, A22, NP, SIZE, flag_ions_first);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end

function [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E)
    % Update stoichiometric submatrix A1
    if temp_NS < temp_NS0
        A11 = eye(temp_NS);
        A11(temp_NG+1:end, temp_NG+1:end) = 0;
        A12 = -A0(temp_ind, temp_ind_E);
        A1 = [A11, A12];
        temp_NS0 = temp_NS;
    end
end

function [temp_ind, temp_ind_swt] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic)
    temp_ind = temp_ind(~ismember(temp_ind, temp_ind_cryogenic));
    temp_ind_swt = temp_ind_swt(~ismember(temp_ind_swt, temp_ind_cryogenic));
end

function A2 = update_matrix_A2(A0_T, A22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E)
    % Update stoichiometric submatrix A2
    A21 = [N0(temp_ind_nswt, 1)' .* A0_T(temp_ind_E, temp_ind_nswt), A0_T(temp_ind_E, temp_ind_swt)];
    A2 = [A21, A22];
end

function A = update_matrix_A(A0_T, A1, A22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E)
    % Update stoichiometric matrix A
    A2 = update_matrix_A2(A0_T, A22, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E);
    A = [A1; A2];
end

function b = update_vector_b(A0, N0, NatomE, ind_E, temp_ind_ions, temp_ind, temp_ind_E, muRT)
    % Update coefficient vector b
    bi = (sum(N0(temp_ind, 1) .* A0(temp_ind, temp_ind_E)))';
    if any(temp_ind_ions)
        bi(temp_ind_E == ind_E) = 0;
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

function [STOP, DeltaN1, Deltab] = compute_STOP(N0, DeltaN0, temp_NG, A0, NatomE, max_NatomE, tolE)
    NPi = sum(N0);
    DeltaN1 = N0 .* abs(DeltaN0) / NPi;
    DeltaN1(temp_NG+1:end) = abs(DeltaN0(temp_NG+1:end)) / NPi;
    Deltab = abs(NatomE - sum(N0 .* A0, 1)');
    Deltab = max(Deltab(NatomE > max_NatomE * tolE));
    STOP = max(max(DeltaN1), Deltab);
end

function [N0, STOP] = check_convergence_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, TOL, TOL_pi, itMax)
    % Check convergence of ionized species
    STOP = 0;
    if any(temp_ind_ions)
        [lambda_ions, ~] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions);
        if abs(lambda_ions) > TOL_pi
            [N0, STOP] = recompute_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, lambda_ions, TOL, TOL_pi, itMax);
        end
    end
end

function [N0, STOP] = recompute_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, lambda_ions, TOL, TOL_pi, itMax)
    % Reestimate composition of ionized species
    
    % Parameters
    A0_ions = A0(temp_ind_ions, ind_E);
    STOP = 1;
    it = 0;
    while STOP > TOL_pi && it < itMax
        it = it + 1;
        % Apply correction
        N0_ions = log(N0(temp_ind_ions, 1)) + A0_ions * lambda_ions;
        % Apply antilog
        N0_ions = exp(N0_ions);
        % Assign values to original vector
        N0(temp_ind_ions, 1) = N0_ions;
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [lambda_ions, ~] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions);
        STOP = abs(lambda_ions);
    end   
    
    Xi_ions = N0(temp_ind_ions, 1) / sum(N0(:, 1));
    if ~any(Xi_ions > TOL)
        STOP = 0;
    end 
end

function [lambda, DeltaN3] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions)
    % Compute relaxation factor ions
    if any(temp_ind_ions)
        lambda = -sum(A0(temp_ind_nswt, ind_E) .* N0(temp_ind_nswt, 1))/ ...
                 sum(A0(temp_ind_nswt, ind_E).^2 .* N0(temp_ind_nswt, 1));
        DeltaN3 = abs(sum(N0(temp_ind_nswt, 1) .* A0(temp_ind_nswt, ind_E)));
    else
        lambda = [];
        DeltaN3 = 0;
    end 
end