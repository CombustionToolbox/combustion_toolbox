function strP = equilibrate_T(self, strR, pP, TP)
    % Set List of Species to List of Products
    self_ListProducts = set_LS_original(self);
    % Compute number of moles 
    [N_ListProducts, DeltaNP, DeltaNP_ions] = select_equilibrium(self_ListProducts, pP, TP, strR);
    % Reshape matrix of number of moles, N
    N = reshape_moles(self, self_ListProducts, N_ListProducts);
    % Compute thermodynamic derivates
    [self.dNi_T, self.dN_T] = equilibrium_dT(self, N, TP, strR);
    [self.dNi_p, self.dN_p] = equilibrium_dp(self, N, strR);
    % Compute properties matrix
    P = SetSpecies(self, self.S.LS, N(:, 1), TP);
    % Compute properties of final mixture
    strP = compute_properties(self, strR, P, pP, TP, DeltaNP, DeltaNP_ions);
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

function pP = compute_pressure(self, strR, TP, N)
    % Compute pressure of product mixture
    pP = (N * TP * self.C.R0 / (strR.v/1e3)) / 1e5;
end

function [N, DeltaNP, DeltaNP_ions] = select_equilibrium(self, pP, TP, strR)
    if ~self.PD.ionization
        % Compute numer of moles without ionization
        [N, DeltaNP] = equilibrium(self, pP, TP, strR);
        DeltaNP_ions = 0;
    else
        % Compute numer of moles with ionization
        [N, DeltaNP, DeltaNP_ions] = equilibrium_ions(self, pP, TP, strR);
    end
end

function strP = compute_properties(self, strR, P, pP, TP, DeltaNP, DeltaNP_ions)
    % Compute properties of final mixture
    if strfind(self.PD.ProblemType, 'P') == 2
        strP = ComputeProperties(self, P, pP, TP);
    else
        NP = sum(P(:, 1) .* (1 - P(:, 10)));
        pP = compute_pressure(self, strR, TP, NP);
        strP = ComputeProperties(self, P, pP, TP);
    end    
    strP.error_moles = DeltaNP;
    strP.error_moles_ions = DeltaNP_ions;
end

function N = reshape_moles(self, self_modified, N_modified)
    % Reshape matrix of number of moles, N
    index = find_ind(self.S.LS, self_modified.S.LS);
    N = self.C.N0.value;
    N(index, 1) = N_modified(:, 1);
end