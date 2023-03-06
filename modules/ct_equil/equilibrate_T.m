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
    %
    % Example:
    %     mix2 = equilibrate(self, self.PS.strR{1}, 1.01325, 3000)
    
    % Check if calculations are for a thermochemical frozen gas (calorically perfect)
    if self.TN.FLAG_TCHEM_FROZEN        
        mix2 = equilibrate_T_tchem(self, mix1, pP, TP);
        return
    end
    
    % Unpack
    guess_moles = unpack(varargin);
    % Check flag
    if ~self.TN.FLAG_FAST, guess_moles = []; end
    % Set List of Species to List of Products
    self_ListProducts = set_LS_original(self, TP);
    % Compute number of moles
    [N, self.dNi_T, self.dN_T, self.dNi_p, self.dN_p, ind_ListProducts, STOP, STOP_ions] = select_equilibrium(self_ListProducts, pP, TP, mix1, guess_moles);
    % Get index of species
    ind = find_ind(self.S.LS, self_ListProducts.S.LS(ind_ListProducts));
    % Reshape composition matrix N, and partial composition partial derivatives 
    N = reshape_vector(self, ind, ind_ListProducts, N);
    self.dNi_T = reshape_vector(self, ind, ind_ListProducts, self.dNi_T);
    self.dNi_p = reshape_vector(self, ind, ind_ListProducts, self.dNi_p);
    % Add moles of frozen species to the moles vector N
    N_mix1 = moles(mix1);
    N(self.S.ind_frozen) = N_mix1(self.S.ind_frozen);
    % Compute property matrix of the species at chemical equilibrium
    % NOTE: If the ind variable is removed from the inputs, the set_species 
    % routine will completely fill the properties matrix
    M0 = set_species(self, self.S.LS, N(:, 1), TP, ind);
    % Compute properties of final mixture
    mix2 = set_properties(self, mix1, M0, pP, TP, STOP, STOP_ions);
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

function self = set_LS_original(self, TP)
    % Set List of Species to List of Products
    
    % Remove ionized species if TP is below T_ions
    if any(self.S.ind_ions) && TP < self.TN.T_ions
        self.Misc.index_LS_original(self.S.ind_ions) = [];
        self.E.ind_E = [];
    end
    % Initialization
    self.S.ind_nswt = []; self.S.ind_swt = [];
    self.S.ind_cryogenic = []; self.S.ind_ions = [];
    % Set list of species for calculations
    self.S.LS = self.S.LS(self.Misc.index_LS_original);
    % Establish cataloged list of species according to the state of the phase
    self = list_phase_species(self, self.S.LS);
    % Update stoichiometric matrix
    self.C.A0.value = self.C.A0.value(self.Misc.index_LS_original, :);
    % Update property matrix
    self.C.M0.value = self.C.M0.value(self.Misc.index_LS_original, :);
    % Update compostion matrix
    self.C.N0.value = self.C.N0.value(self.Misc.index_LS_original, :);
end

function pP = set_pressure(self, mix1, TP, N)
    % Compute pressure of product mixture
    if strcmpi(self.PD.ProblemType, 'SV')
        pP = mix1.p * self.PD.EOS.pressure(self, N / mix1.N, TP / mix1.T, self.PD.vP_vR.value, self.S.LS, mix1.Xi);
    else
        pP = self.PD.EOS.pressure(self, N, TP, mix1.v, self.S.LS, mix1.Xi) * 1e-5;
    end

end

function [N, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions] = select_equilibrium(self, pP, TP, mix1, guess_moles)
    % Select equilibrium: TP: Gibbs; TV: Helmholtz
    if strfind(self.PD.ProblemType, 'P') == 2
        [N, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions] = equilibrium_gibbs(self, pP, TP, mix1, guess_moles);
    else
        [N, dNi_T, dN_T, dNi_p, dN_p, ind, STOP, STOP_ions] = equilibrium_helmholtz(self, mix1.v, TP, mix1, guess_moles);
    end

end

function mix2 = set_properties(self, mix1, properties_matrix, pP, TP, STOP, STOP_ions)
    % Compute properties of final mixture
    if strfind(self.PD.ProblemType, 'P') == 2
        mix2 = compute_properties(self, properties_matrix, pP, TP);
    else
        NP = sum(properties_matrix(:, self.C.M0.ind_ni) .* (1 - properties_matrix(:, self.C.M0.ind_phase)));
        pP = set_pressure(self, mix1, TP, NP);
        mix2 = compute_properties(self, properties_matrix, pP, TP);
    end

    mix2.error_moles = STOP;
    mix2.error_moles_ions = STOP_ions;
end

function vector = reshape_vector(self, index, ind_modified, vector_modified)
    % Reshape vector containing all the species
    vector = self.C.N0.value(:, 1);
    vector(index, 1) = vector_modified(ind_modified, 1);
end