function self = check_inputs(self)
    % Check that all the inputs are specified
    if ~self.Misc.FLAG_CHECK_INPUTS
        check_inputs_prop(self, 'TR');
        check_inputs_prop(self, 'pR');
        check_inputs_species(self);
        switch self.PD.ProblemType
            case 'TP' % * TP: Equilibrium composition at defined T and p
                check_inputs_prop(self, 'TP');
                check_inputs_prop(self, 'pP');
            case 'HP' % * HP: Adiabatic T and composition at constant p
                check_inputs_prop(self, 'pP');
            case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
                check_inputs_prop(self, 'pP');
            case 'TV' % * TV: Equilibrium composition at defined T and constant v
                check_inputs_prop(self, 'TP');
                self = set_prop(self, 'pP', self.PD.pR.value); % Guess
            case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
                self = set_prop(self, 'pP', self.PD.pR.value); % Guess
            case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
                check_inputs_prop(self, 'vP_vR');
            case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
                check_inputs_prop(self, 'u1');
            case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
                check_inputs_prop(self, 'u1');
            case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION
                check_inputs_prop(self, 'overdriven');
        end
    end
    self.Misc.FLAG_CHECK_INPUTS = true;
end

% SUB-PASS FUNCTIONS
function check_inputs_prop(self, name)
    % Check that the property is specified
    if isempty(self.PD.(name).value)
        error('ERROR: %s is not specified', self.PD.(name).description);
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