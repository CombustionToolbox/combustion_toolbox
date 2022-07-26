function mix2 = equilibrate_T(self, mix1, pP, TP, varargin)
    % Obtain equilibrium properties and composition for the given temperature [K] and pressure [bar]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float): Pressure [bar]
    %     TP (float): Temperature [K]
    %
    % Optional Args:
    %     guess_moles (float): mixture composition [mol] of a previous computation
    %     
    % Returns:
    %     mix2 (struct): Properties of the final mixture
    
    % Definitions
    FLAG_FAST = self.TN.FLAG_FAST;
    % Unpack
    guess_moles = unpack(varargin);
    % Check flag
    if ~FLAG_FAST, guess_moles = []; end
    % Set List of Species to List of Products
    self_ListProducts = set_LS_original(self);
    % Compute number of moles 
    [N, self.dNi_T, self.dN_T, self.dNi_p, self.dN_p, STOP, STOP_ions] = select_equilibrium(self_ListProducts, pP, TP, mix1, guess_moles);
    % Reshape matrix of number of moles, N
    self.dNi_T = reshape_vectors(self, self_ListProducts, self.dNi_T);
    self.dNi_p = reshape_vectors(self, self_ListProducts, self.dNi_p);
    % Add moles of frozen species to the moles vector N
    N_mix1 = moles(mix1);
    N(self.S.ind_frozen) = N_mix1(self.S.ind_frozen);
    % Compute properties matrix
    P = SetSpecies(self, self.S.LS, N(:, 1), TP);
    % Compute properties of final mixture
    mix2 = compute_properties(self, mix1, P, pP, TP, STOP, STOP_ions);
end

% SUB-PASS FUNCTIONS
function guess_moles = unpack(value)
    % Unpack inputs
    if isempty(value)
        guess_moles = [];
    else
        guess_moles = value{1};
    end
end

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
        pP = mix1.p * self.PD.EOS.pressure(self, N/mix1.N, TP/mix1.T, self.PD.vP_vR.value);
    else
        pP = self.PD.EOS.pressure(self, N, TP, mix1.v) * 1e-5;
    end
end

function [N, dNi_T, dN_T, dNi_p, dN_p, STOP, STOP_ions] = select_equilibrium(self, pP, TP, mix1, guess_moles)
    % Select equilibrium: TP: Gibbs; TV: Helmholtz
    if strfind(self.PD.ProblemType, 'P') == 2
        [N, dNi_T, dN_T, dNi_p, dN_p, STOP, STOP_ions] = equilibrium(self, pP, TP, mix1, guess_moles);
    else
        [N, dNi_T, dN_T, dNi_p, dN_p, STOP, STOP_ions] = equilibrium_helmholtz(self, mix1.v, TP, mix1, guess_moles);
    end
end

function mix2 = compute_properties(self, mix1, P, pP, TP, STOP, STOP_ions)
    % Compute properties of final mixture
    if strfind(self.PD.ProblemType, 'P') == 2
        mix2 = ComputeProperties(self, P, pP, TP);
    else
        NP = sum(P(:, self.C.M0.ind_ni) .* (1 - P(:, self.C.M0.ind_phase)));
        pP = compute_pressure(self, mix1, TP, NP);
        mix2 = ComputeProperties(self, P, pP, TP);
    end    
    mix2.error_moles = STOP;
    mix2.error_moles_ions = STOP_ions;
end

function vector = reshape_vectors(self, self_modified, vector_modified)
    % Reshape vector containing all the species
    index = find_ind(self.S.LS, self_modified.S.LS);
    vector = self.C.N0.value(:, 1);
    vector(index, 1) = vector_modified(:, 1);
end