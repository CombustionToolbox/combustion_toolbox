function self = check_inputs(self)
    % Check that all the inputs are specified
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases
    
    if self.Misc.FLAG_CHECK_INPUTS
        return
    end

    self.C.l_phi = length(self.PD.phi.value);

    if isempty(self.PD.TR.value)
        self = set_prop(self, 'TR', 300); % In case is not specified (species with different temperature)
    else
        self = check_inputs_prop(self, 'TR');
    end

    self = check_inputs_prop(self, 'pR');
    self = check_inputs_species(self);

    switch self.PD.ProblemType
        case 'TP' % * TP: Equilibrium composition at defined T and p
            self = check_inputs_prop(self, 'TP');
            self = check_inputs_prop(self, 'pP');
        case 'HP' % * HP: Adiabatic T and composition at constant p
            self = check_inputs_prop(self, 'pP');
        case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            self = check_inputs_prop(self, 'pP');
        case 'TV' % * TV: Equilibrium composition at defined T and constant v
            self = check_inputs_prop(self, 'TP');
            self = set_prop(self, 'pP', self.PD.pR.value); % Guess
        case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
            self = set_prop(self, 'pP', self.PD.pR.value); % Guess
        case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            self = check_inputs_prop(self, 'vP_vR');
            self = set_prop(self, 'pP', self.PD.pR.value); % Guess
        case {'SHOCK_I', 'SHOCK_R'} % * SHOCK_I and SHOCK_R: Calculate planar shock wave
            
            % Check pre-shock velocity/Mach
            try
                self = check_inputs_prop(self, 'u1');
            catch
                self = set_velocity_from_Mach(self);
            end

        case {'SHOCK_POLAR', 'SHOCK_POLAR_LIMITRR'} % * SHOCK_POLAR: Calculate oblique shock polars
            
            % Check pre-shock velocity/Mach
            try
                self = check_inputs_prop(self, 'u1');
            catch
                self = set_velocity_from_Mach(self);
            end

        case {'SHOCK_OBLIQUE', 'SHOCK_POLAR_R'} % * SHOCK_OBLIQUE OR SHOCK_POLAR
            
            % Check pre-shock velocity/Mach
            try
                self = check_inputs_prop(self, 'u1');
            catch
                self = set_velocity_from_Mach(self);
            end
            
            try
                % Compute from deflection angle (1 solution)
                self = check_inputs_prop(self, 'theta'); 
            catch
                % Compute from wave angle (2 solutions: weak and stron shocks)
                self = check_inputs_prop(self, 'beta');
            end

        case {'DET_OBLIQUE'} % * DET_OBLIQUE: Compute oblique detonation
            self = check_inputs_prop(self, 'drive_factor');

            try
                % Compute from deflection angle (1 solution)
                self = check_inputs_prop(self, 'theta');
            catch
                % Compute from wave angle (2 solutions: weak and stron shocks)
                self = check_inputs_prop(self, 'beta');
            end

        case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R', 'DET_POLAR'} % * DET_OVERDRIVEN, DET_UNDERDRIVEN, DET_POLAR, incident and reflected states
            self = check_inputs_prop(self, 'drive_factor');
        case {'ROCKET'} % * ROCKET: Rocket propellant performance

            if ~self.PD.FLAG_IAC
                self = check_inputs_prop(self, 'Aratio_c');
            end

            if ~isempty(self.PD.Aratio.value)
                self = check_inputs_prop(self, 'Aratio');
            end

    end

    % Check if a stoichiometric equivalence ratio is considered
    if any(self.PD.phi.value == 1)
        self.PD.phi.value(self.PD.phi.value == 1) = 1 + self.TN.tolN * 10;
    end

    % Check reactant species are contained in the list of products (initial computations)
    self = check_FOI(self, [self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert]);

    % Update FLAG
    self.Misc.FLAG_CHECK_INPUTS = true;
end

% SUB-PASS FUNCTIONS
function self = check_inputs_prop(self, name)
    % Check that the property is specified
    if isempty(self.PD.(name).value)
        error('ERROR: %s is not specified', self.PD.(name).description);
    end

    value_length = self.PD.(name).value;
    self = set_length_phi(self, name, value_length);
end

function self = set_length_phi(self, name, value)
    % Set length equivalence ratio
    if self.Misc.FLAG_LENGTH_LOOP
        return
    end

    if self.Misc.FLAGS_PROP.phi
        self.C.l_phi = length(self.PD.phi.value);
        self.PD.range = self.PD.phi.value;
        self.PD.range_name = 'phi';
        self.Misc.FLAG_LENGTH_LOOP = true;
        return
    end

    if self.Misc.FLAGS_PROP.(name)
        
        if isempty(self.PD.phi.value)
            self.PD.phi.value = 1;
        end

        self.PD.phi.value = self.PD.phi.value(1) * ones(length(value));
        self.Misc.FLAG_LENGTH_LOOP = true;
    elseif ~isempty(self.PD.phi.value)
        self.PD.phi.value = self.PD.phi.value(1);
    else
        self.PD.phi.value = 1;
    end

    self.C.l_phi = length(self.PD.phi.value);
    self.PD.range = value;
    self.PD.range_name = name;
end

function self = check_inputs_species(self)
    % Check that species, the Nº moles, and the temperature of the species
    % are specified
    if isempty(self.PD.S_Fuel) && isempty(self.PD.S_Oxidizer) && isempty(self.PD.S_Inert)
        error('ERROR: species are not specified');
    end

    if isempty(self.PD.N_Fuel) && isempty(self.PD.N_Oxidizer) && isempty(self.PD.N_Inert)
        fprintf('COMPUTING Nº MOLES FROM EQUIVALENCE RATIO (PHI).\n');
    end

    if length(self.PD.S_Fuel) ~= length(self.PD.S_Fuel)
        error('ERROR: mismatch length fuel species and Nº moles');
    end

    if length(self.PD.S_Oxidizer) ~= length(self.PD.S_Oxidizer)
        error('ERROR: mismatch length oxidizer species and Nº moles');
    end

    if length(self.PD.S_Inert) ~= length(self.PD.S_Inert)
        error('ERROR: mismatch length inert species and Nº moles');
    end

    if length(self.PD.TR.value) == 1
        self = check_input_temperatures(self, 'Fuel');
        self = check_input_temperatures(self, 'Oxidizer');
        self = check_input_temperatures(self, 'Inert');
    end

end

function self = check_input_temperatures(self, name)
    % Assign temperature of the species in case it was not defined
    Sname = strcat('S_', name);
    Tname = strcat('T_', name);

    if isempty(self.PD.(Sname))
        return
    end
    
    if ~isempty(self.PD.(Tname))
        return
    end

    species = self.PD.(Sname);
    N = length(species);

    for i = N:-1:1

        if self.DB.(species{i}).phase

            if self.PD.TR.value >= min(self.DB.(species{i}).T) && self.PD.TR.value <= max(self.DB.(species{i}).T)
                T(i) = self.PD.TR.value;
            else
                T(i) = self.DB.(species{i}).T;
            end

        else
            T(i) = self.PD.TR.value;
        end

    end

    self.PD.(Tname) = T;
end

function self = set_velocity_from_Mach(self)
    % Set pre-shock velocity from pre-shock Mach number
    % 
    % Note:
    %     This is only valid with fixed reactant's composition, because
    %     check_inputs routines only perform this check one time

    self = check_inputs_prop(self, 'M1');
    species = [self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert];
    moles = [self.PD.N_Fuel, self.PD.N_Oxidizer, self.PD.N_Inert];
    sound_velocity = compute_sound(self.PD.TR.value, self.PD.pR.value, species, moles, 'self', self);
    u1 = self.PD.M1.value * sound_velocity;
    self = set_prop(self, 'u1', u1);
end