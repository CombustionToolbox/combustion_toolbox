function self = Constants()
    % Initialize struct with constants data
    % 
    % Returns:
    %     self (struct): struct with constants data

    % Description
    self.description = "Constants";
    % Variables
    [self.release, self.date] = get_combustion_toolbox_version();
    self.R0 = 8.31446261815324; % [J/(K mol)]. Universal gas constant
    self.gravity = 9.80665;     % [m/s2]. Standard gravity
    self.A0.description = "Stoichiometric Matrix: number of atoms of each element contained in each species";
    self.A0.value = [];
    self.M0.description = "Matrix with properties of each species";
    self.M0.value = [];
    self.M0.ind_hfi   = 1; % Index enthalpy of formation
    self.M0.ind_efi   = 2; % Index internal energy of formation
    self.M0.ind_W     = 3; % Index molecular weight
    self.M0.ind_phase = 4; % Index phase
    self.M0.ind_ni    = 5; % Index number of moles
    self.M0.ind_hi    = 6; % Index enthalpy
    self.M0.ind_cPi   = 7; % Index specific heat at constant pressure
    self.M0.ind_si    = 8; % Index entropy
    self.N_prop.description = 'Number of properties in properties_matrix';
    self.N_prop.value = 8;
    self.N0.description = "Reduced Matrix with number of moles and phase of each species";
    self.N0.value = [];
    self.MassorMolar = 'mass';
    self.firstrow = true;
    self.mintol_display = 1e-14;
    self.l_phi = []; % length phi vector
    self.composition_units = 'molar fraction'; % Possible values: mol, molar fraction or mass fraction
end
