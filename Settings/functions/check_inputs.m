function self = check_inputs(self)
    % Check that all the inputs are specified
    if ~self.Misc.FLAG_CHECK_INPUTS
        self.C.l_phi = length(self.PD.phi.value);
        self = check_inputs_prop(self, 'TR');
        self = check_inputs_prop(self, 'pR');
        check_inputs_species(self);
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
            case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
                self = check_inputs_prop(self, 'u1');
            case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
                self = check_inputs_prop(self, 'u1');
            case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION
                self = check_inputs_prop(self, 'overdriven');
        end
    end
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
    if ~self.Misc.FLAG_LENGTH_LOOP
        if isfield(self.Misc.FLAGS_PROP, 'phi')
            if self.Misc.FLAGS_PROP.phi
                self.C.l_phi = length(self.PD.phi.value);
                self.PD.range = self.PD.phi.value;
                self.PD.range_name = 'phi';
                self.Misc.FLAG_LENGTH_LOOP = true;
                return
            end
        end
        if self.Misc.FLAGS_PROP.(name)
            self.PD.phi.value = self.PD.phi.value(1) * ones(length(value));
            self.Misc.FLAG_LENGTH_LOOP = true;
        else
            self.PD.phi.value = 1;
        end
        self.C.l_phi = length(self.PD.phi.value);
        self.PD.range = value;
        self.PD.range_name = name;
    end
end

function check_inputs_species(self) 
    % Check that species and the Nº moles are specified
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
end