function mix =  ComputeProperties(self, SpeciesMatrix, p, T)
    % Compute properties from the given SpeciesMatrix at pressure p [bar] 
    % and temperature T [K]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     SpeciesMatrix (float): Matrix with the properties values of the mixture
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
    % Unpack SpeciesMatrix
    Ni          = SpeciesMatrix(:,1);                       % [mol]
    mix.N       = sum(SpeciesMatrix(:, self.C.M0.ind_ni));  % [mol] 
    mix.hf      = sum(SpeciesMatrix(:, self.C.M0.ind_hfi)); % [kJ]
    mix.h       = sum(SpeciesMatrix(:, self.C.M0.ind_hi));  % [kJ]
    mix.ef      = sum(SpeciesMatrix(:, self.C.M0.ind_efi)); % [kJ]
    mix.cP      = sum(SpeciesMatrix(:, self.C.M0.ind_cPi)); % [J/K]
    mix.S0      = sum(SpeciesMatrix(:, self.C.M0.ind_si));  % [kJ/K]
    mix.pv      = sum(SpeciesMatrix(:, self.C.M0.ind_pVi)); % [bar m3]
    mix.phase   = SpeciesMatrix(:, self.C.M0.ind_phase);    % [bool]
    mix.mi      = sum(SpeciesMatrix(:, self.C.M0.ind_mi));  % [kg]
    % Compute mass fractions [-]
    mix.Yi = SpeciesMatrix(:, self.C.M0.ind_mi)./mix.mi;
    % Compute mean molecular weight [g/mol]
    mix.W = 1/sum(mix.Yi./SpeciesMatrix(:, self.C.M0.ind_mmi), 'OmitNan');
    % Compute molar fractions [-]
    mix.Xi = Ni/mix.N;
    % Get non zero species
    FLAG_NONZERO = mix.Xi > 0;
    % Compute vector atoms of each element
    mix.NatomE = sum(Ni .* self.C.A0.value);
    % Compute vector atoms of each element without frozen species
    mix.NatomE_react = sum(SpeciesMatrix(self.S.ind_react, self.C.M0.ind_ni) .* self.C.A0.value(self.S.ind_react, :));
    % Compute volume [m3]
    mix.v = mix.pv / mix.p;
    % Compute density [kg/m3]
    mix.rho = mix.mi / mix.v;
    % Compute internal energy [kJ]
    mix.e = mix.h - molesGas(mix) * R0 * T * 1e-3;
    % Compute thermal internal energy [kJ]
    mix.DeT = mix.e - mix.ef;
    % Compute thermal enthalpy [kJ]
    mix.DhT = mix.h - mix.hf;
    % Compute entropy of mixing [kJ/K]
    mix.DS = compute_entropy_mixing(mix, Ni, R0, FLAG_NONZERO);
    % Compute entropy [kJ/K]
    mix.S = mix.S0 + mix.DS;
    % Compute Gibbs energy [kJ]
    mix.g = mix.h - mix.T * mix.S;
    % Compute specific heat at constant volume [J/K]
    mix.cV = mix.cP - R0 * mix.N;
    % Compute Adibatic index [-]
    mix.gamma = mix.cP/mix.cV;
    % Compute sound velocity [m/s]
    mix.sound = sqrt(mix.gamma * p * 1e5 / mix.rho);
    
    % Correction of: cP, cV, gamma, and speed of sound as consequence of the
    % chemical reaction
    if isfield(self, 'dNi_T')
        mix.dVdT_p =  1 + self.dN_T; % [-]
        mix.dVdp_T = -1 + self.dN_p; % [-]
        if ~any(isnan(self.dNi_T)) && ~any(isinf(self.dNi_T))
            delta = ~mix.phase;
            h0_j = (SpeciesMatrix(:, self.C.M0.ind_hi)) ./ Ni * 1e3; % [J/mol]
            mix.cP_r = sum(h0_j/T .* (1 +  delta .* (Ni - 1)) .* self.dNi_T, 'omitnan'); % [J/K]
            mix.cP_f = mix.cP;
            mix.cP = mix.cP_f + mix.cP_r; % [J/K]
            mix.cV = mix.cP + (molesGas(mix) * R0 * mix.dVdT_p^2) / mix.dVdp_T; % [J/K]
            mix.gamma = mix.cP/mix.cV; % [-]
            mix.gamma_s = - mix.gamma / mix.dVdp_T; % [-]
            mix.sound = sqrt(mix.gamma_s * p*1e5 / mix.rho); % [m/s]
        else
            mix.gamma_s = - 1 / mix.dVdp_T; % [-]
        end
    else
        mix.gamma_s = mix.gamma;
    end
end

% SUB-PASS FUNCTIONS
function DS = compute_entropy_mixing(mix, Ni, R0, FLAG_NONZERO)
    % Compute entropy of mixing [kJ/K]. 
    % Note: only nonzero for gaseous species

    DSi = Ni(FLAG_NONZERO) .* log(Ni(FLAG_NONZERO) / molesGas(mix) * mix.p) .* (1 - mix.phase(FLAG_NONZERO));
    DS  = -R0 * sum(DSi) * 1e-3;
end