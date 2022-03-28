function mix2 = equilibrate_T(self, mix1, pP, TP)
    % Obtain equilibrium properties and composition for the given temperature [K] and pressure [bar]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float):    Pressure [bar]
    %     TP (float):    Temperature [K]
    %
    % Returns:
    %     mix2 (struct): Properties of the final mixture

    % Set List of Species to List of Products
    self_ListProducts = set_LS_original(self);
    % Compute number of moles 
    [N_ListProducts, DeltaNP, DeltaNP_ions] = select_equilibrium(self_ListProducts, pP, TP, mix1);
    % Reshape matrix of number of moles, N
    N = reshape_moles(self, self_ListProducts, N_ListProducts);
    % Compute thermodynamic derivates
    [self.dNi_T, self.dN_T] = equilibrium_dT(self, N, TP, mix1);
    [self.dNi_p, self.dN_p] = equilibrium_dp(self, N, mix1);
    % Compute properties matrix
    P = SetSpecies(self, self.S.LS, N(:, 1), TP);
    % Compute properties of final mixture
    mix2 = compute_properties(self, mix1, P, pP, TP, DeltaNP, DeltaNP_ions);
end

% SUB-PASS FUNCTIONS
function self = set_LS_original(self)
    % Set List of Species to List of Products
    self.S.ind_nswt = []; self.S.ind_swt = []; self.S.ind_cryogenic = [];
    self.S.LS = self.S.LS(self.Misc.index_LS_original);
    self = list_phase_species(self, self.S.LS);
    self.C.A0.value = self.C.A0.value(self.Misc.index_LS_original, :);
    self.C.M0.value = self.C.M0.value(self.Misc.index_LS_original, :);
    self.C.N0.value = self.C.N0.value(self.Misc.index_LS_original, :);
end

function pP = compute_pressure(self, mix1, TP, N)
    % Compute pressure of product mixture
    if strcmpi(self.PD.ProblemType, 'SV')
        pP = mix1.p * N/mix1.N * TP/mix1.T * self.PD.vP_vR.value;
    else
        pP = (N * TP * self.C.R0 / mix1.v) * 1e-5;
    end
end

function [N, DeltaNP, DeltaNP_ions] = select_equilibrium(self, pP, TP, mix1)
    if ~self.PD.ionization
        % Compute numer of moles without ionization
        [N, DeltaNP] = equilibrium(self, pP, TP, mix1);
        DeltaNP_ions = 0;
    else
        % Compute numer of moles with ionization
        [N, DeltaNP, DeltaNP_ions] = equilibrium_ions(self, pP, TP, mix1);
    end
end

function mix2 = compute_properties(self, mix1, P, pP, TP, DeltaNP, DeltaNP_ions)
    % Compute properties of final mixture
    if strfind(self.PD.ProblemType, 'P') == 2
        mix2 = ComputeProperties(self, P, pP, TP);
    else
        NP = sum(P(:, 1) .* (1 - P(:, 10)));
        pP = compute_pressure(self, mix1, TP, NP);
        mix2 = ComputeProperties(self, P, pP, TP);
    end    
    mix2.error_moles = DeltaNP;
    mix2.error_moles_ions = DeltaNP_ions;
end

function N = reshape_moles(self, self_modified, N_modified)
    % Reshape matrix of number of moles, N
    index = find_ind(self.S.LS, self_modified.S.LS);
    N = self.C.N0.value;
    N(index, 1) = N_modified(:, 1);
end