function mix = compute_properties(self, properties_matrix, p, T)
    % Compute properties from the given properties matrix at pressure p [bar]
    % and temperature T [K]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     properties_matrix (float): Matrix with the properties of the mixture
    %     p (float): Pressure [bar]
    %     T (float): Temperature [K]
    %
    % Returns:
    %     mix (struct): Properties of the mixture

    % Definitions
    R0 = self.C.R0; % [J/(K mol)] Universal gas constant
    % Inputs
    mix.p = p; % [bar]
    mix.T = T; % [K]
    % Unpack properties_matrix
    Ni = properties_matrix(:, self.C.M0.ind_ni); % [mol]
    mix.N = sum(properties_matrix(:, self.C.M0.ind_ni)); % [mol]
    mix.hf = dot(properties_matrix(:, self.C.M0.ind_hfi), Ni); % [kJ]
    mix.h = dot(properties_matrix(:, self.C.M0.ind_hi), Ni); % [kJ]
    mix.ef = dot(properties_matrix(:, self.C.M0.ind_efi), Ni); % [kJ]
    mix.cP = dot(properties_matrix(:, self.C.M0.ind_cPi), Ni); % [J/K]
    mix.S0 = dot(properties_matrix(:, self.C.M0.ind_si), Ni); % [kJ/K]
    mix.phase = properties_matrix(:, self.C.M0.ind_phase); % [bool]
    % Compute total composition of gas species [mol]
    N_gas = sum(Ni(mix.phase == 0));
    % Compute molar fractions [-]
    mix.Xi = Ni / mix.N;
    % Compute molecular weight [g/mol]
    mix.W = dot(Ni, properties_matrix(:, self.C.M0.ind_W)) / N_gas;
    % Compute mean molecular weight [g/mol]
    mix.MW = dot(Ni, properties_matrix(:, self.C.M0.ind_W)) / mix.N;
    % Compute mass mixture [kg]
    mix.mi = mix.MW * mix.N * 1e-3; % [kg]
    % Compute mass fractions [-]
    mix.Yi = (Ni .* properties_matrix(:, self.C.M0.ind_W) * 1e-3) ./ mix.mi;
    % Get non zero species
    FLAG_NONZERO = mix.Xi > 0;
    % Compute vector atoms of each element
    mix.NatomE = sum(Ni .* self.C.A0.value, 1);
    % Compute vector atoms of each element without frozen species
    mix.NatomE_react = sum(properties_matrix(self.S.ind_react, self.C.M0.ind_ni) .* self.C.A0.value(self.S.ind_react, :), 1);
    % Compute volume [m3]
    mix.v = self.PD.EOS.volume(self, T, convert_bar_to_Pa(p), self.S.LS, mix.Xi) * N_gas;
    % Compute density [kg/m3]
    mix.rho = mix.mi / mix.v;
    % Compute internal energy [kJ]
    mix.e = mix.h - N_gas * R0 * T * 1e-3;
    % Compute thermal internal energy [kJ]
    mix.DeT = mix.e - mix.ef;
    % Compute thermal enthalpy [kJ]
    mix.DhT = mix.h - mix.hf;
    % Compute entropy of mixing [kJ/K]
    mix.DS = compute_entropy_mixing(mix, Ni, N_gas, R0, FLAG_NONZERO);
    % Compute entropy [kJ/K]
    mix.S = mix.S0 + mix.DS;
    % Compute Gibbs energy [kJ]
    mix.g = mix.h - mix.T * mix.S;
    % Compute specific heat at constant volume [J/K]
    mix.cV = mix.cP - R0 * N_gas;
    % Compute Adibatic index [-]
    mix.gamma = mix.cP / mix.cV;
    % Compute sound velocity [m/s]
    mix.sound = sqrt(mix.gamma * convert_bar_to_Pa(p) / mix.rho);

    % Correction of: cP, cV, gamma, and speed of sound as consequence of the
    % chemical reaction
    if isfield(self, 'dNi_T')
        mix.dVdT_p = 1 + self.dN_T; % [-]
        mix.dVdp_T = -1 + self.dN_p; % [-]

        if ~any(isnan(self.dNi_T)) && ~any(isinf(self.dNi_T))
            delta = ~mix.phase;
            h0_j = properties_matrix(:, self.C.M0.ind_hi) * 1e3; % [J/mol]
            mix.cP_r = sum(h0_j / T .* (1 + delta .* (Ni - 1)) .* self.dNi_T, 'omitnan'); % [J/K]
            mix.cP_f = mix.cP;
            mix.cP = mix.cP_f + mix.cP_r; % [J/K]
            mix.cV = mix.cP + (N_gas * R0 * mix.dVdT_p^2) / mix.dVdp_T; % [J/K]
            mix.gamma = mix.cP / mix.cV; % [-]
            mix.gamma_s =- mix.gamma / mix.dVdp_T; % [-]
            mix.sound = sqrt(mix.gamma_s * convert_bar_to_Pa(p) / mix.rho); % [m/s]
        else
            mix.gamma_s =- 1 / mix.dVdp_T; % [-]
        end

    else
        mix.gamma_s = mix.gamma;
    end

end

% SUB-PASS FUNCTIONS
function DS = compute_entropy_mixing(mix, Ni, N_gas, R0, FLAG_NONZERO)
    % Compute entropy of mixing [kJ/K].
    % Note: only nonzero for gaseous species
    DSi = Ni(FLAG_NONZERO) .* log(Ni(FLAG_NONZERO) / N_gas * mix.p) .* (1 - mix.phase(FLAG_NONZERO));
    DS = -R0 * sum(DSi) * 1e-3;
end
