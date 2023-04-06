function self = Constants()
    % Initialize struct with constants data
    %
    % Attributes:
    %     description (char): Description of the struct
    %     release (char): Release of the Combustion Toolbox
    %     date (char): Date of the release
    %     R0 (float): Universal gas constant [J/(K mol)]
    %     gravity (float): Standard gravity [m/s2]
    %     A0 (struct): Stoichiometric Matrix
    %     M0 (struct): Matrix with properties of each species
    %     N_prop (struct): Number of properties in properties_matrix
    %     N0 (struct): Reduced Matrix with number of moles and phase of each species
    %     MassorMolar (char): 'mass' or 'molar'
    %     mintol_display (float): Minimum tolerance to display results
    %     l_phi (float): Length equivalence ratio vector
    %     composition_units (char): Possible values: mol, molar fraction or mass fraction
    %
    % Returns:
    %     self (struct): Struct with constants data

    % Description
    self.description = 'Constants';
    % Variables
    [self.release, self.date] = get_combustion_toolbox_version();
    self.R0 = 8.31446261815324; % [J/(K mol)]. Universal gas constant
    self.gravity = 9.80665;     % [m/s2]. Standard gravity
    self.A0.description = 'Stoichiometric Matrix: number of atoms of each element contained in each species';
    self.A0.value = [];
    self.M0.description = 'Matrix with properties of each species';
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
    self.N0.description = 'Reduced Matrix with number of moles and phase of each species';
    self.N0.value = [];
    self.MassorMolar = 'mass';
    self.firstrow = true; % (will be deprecated)
    self.mintol_display = 1e-14; % (will be moved to Miscellaneous)
    self.l_phi = []; % length equivalence ratio vector (will be deprecated)
    self.composition_units = 'molar fraction'; % Possible values: mol, molar fraction or mass fraction
end
