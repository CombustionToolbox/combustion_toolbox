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
    Ni          = SpeciesMatrix(:,1);       % [mol]
    mix.N       = sum(SpeciesMatrix(:,1));  % [mol] 
    mix.hf      = sum(SpeciesMatrix(:,2));  % [kJ]
    mix.DhT     = sum(SpeciesMatrix(:,3));  % [kJ]
    mix.ef      = sum(SpeciesMatrix(:,4));  % [kJ]
    mix.DeT     = sum(SpeciesMatrix(:,5));  % [kJ]
    mix.cP      = sum(SpeciesMatrix(:,6));  % [J/K]     
    mix.cV      = sum(SpeciesMatrix(:,7));  % [J/K]  
    mix.S0      = sum(SpeciesMatrix(:,8));  % [kJ/K]
    mix.pv      = sum(SpeciesMatrix(:,9));  % [bar m3]
    mix.swtCond = SpeciesMatrix(:,10);      % [bool]
    mix.mi      = sum(SpeciesMatrix(:,11)); % [kg]
    % Compute mass fractions [-]
    mix.Yi = SpeciesMatrix(:,11)./mix.mi; % [-]
    % Compute molar fractions [-]
    mix.Xi = Ni/mix.N;
    % Get non zero species
    FLAG_NONZERO = mix.Xi > 0; 
    % Compute mean molecular weight [g/mol]
    mix.W = 1/sum(mix.Yi./SpeciesMatrix(:,12), 'OmitNan');
    % Compute vector atoms of each element
    mix.NatomE = sum(Ni .* self.C.A0.value);
    % Assign values for C, H, O, and N elements
    mix = assign_values_CHON(self, mix);
    % Compute volume [m3]
    mix.v = mix.pv / mix.p;
    % Compute density [kg/m3]
    mix.rho = mix.mi / mix.v;
    % Compute enthalpy [kJ]
    mix.h = mix.hf + mix.DhT;
    % Compute internal energy [kJ]
    mix.e = mix.ef + mix.DeT; 
    % Compute entropy of mixing [kJ/K]
    mix.DS = compute_entropy_mixing(mix, Ni, R0, FLAG_NONZERO);
    % Compute entropy [kJ/K]
    mix.S = mix.S0 + mix.DS;
    % Compute Gibbs energy [kJ]
    mix.g = mix.h - mix.T * mix.S;
    % Compute internal energy [kJ]
    mix.e = mix.h - sum(Ni(FLAG_NONZERO) .* (1 - mix.swtCond(FLAG_NONZERO))) * R0 * T *1e-3; % [kJ]
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
            delta = ~mix.swtCond;
            h0_j = (SpeciesMatrix(:, 2) + SpeciesMatrix(:, 3)) ./ Ni * 1e3; % [J/mol]
            mix.cP_r = sum(h0_j/T .* (1 +  delta .* (Ni - 1)) .* self.dNi_T, 'omitnan'); % [J/K]
            mix.cP = mix.cP + mix.cP_r; % [J/K]
            mix.cV = mix.cP + (mix.pv/T * mix.dVdT_p^2) / mix.dVdp_T * 1e5; % [J/K]
            mix.gamma = mix.cP/mix.cV; % [-]
            mix.gamma_s = - mix.gamma / mix.dVdp_T; % [-]
            mix.sound = sqrt(mix.gamma_s*p*1e5/mix.rho); % [m/s]
        else
            mix.gamma_s = - 1 / mix.dVdp_T; % [-]
        end
    else
        mix.gamma_s = mix.gamma;
    end
end

% SUB-PASS FUNCTIONS
function mix = assign_values_CHON(self, mix)
    % Assign values for C, H, O, and N elements
    
    if isempty(self.E.ind_C), mix.x = 0; else, mix.x = mix.NatomE(self.E.ind_C); end
    if isempty(self.E.ind_H), mix.y = 0; else, mix.y = mix.NatomE(self.E.ind_H); end
    if isempty(self.E.ind_O), mix.z = 0; else, mix.z = mix.NatomE(self.E.ind_O); end
    if isempty(self.E.ind_N), mix.w = 0; else, mix.w = mix.NatomE(self.E.ind_N); end
end

function DS = compute_entropy_mixing(mix, Ni, R0, FLAG_NONZERO)
    % Compute entropy of mixing [kJ/K].
    % Note: only nonzero for noncondensed species

    DSi = Ni(FLAG_NONZERO) .* log(mix.Xi(FLAG_NONZERO) * mix.p) .* (1 - mix.swtCond(FLAG_NONZERO));
    DS  = -R0 * sum(DSi) * 1e-3;
end