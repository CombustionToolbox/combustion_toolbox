function gui_ReactantsValueChanged(self)
    % Update reactants & GUI
    self = gui_update_Reactants(self);
    % Update UITable_R2 (species, numer of moles, mole fractions)
    self.UITable_R2.Data = self.UITable_R.Data(:, 1:3); 
    % Get temperature of the mixture
    self.UITable_R.Data(:, 5) = {sscanf(self.PR1.Value, '%f')};
    % Update equivalence ratio 
    self.edit_phi2.Value = self.edit_phi.Value;
end

% SUB-PASS FUNCTIONS
function self = gui_update_Reactants(self)
    switch self.Reactants.Value
        case '1' % No species selected
            self.UITable_R.Data  = [];
            self.UITable_P.Data  = [];
            self.UITable_R2.Data = [];
            return
        case '2' % AIR IDEAL (21% O2 + 79% N2)
            self.UITable_R.Data = {'N2', 0.79/0.21, 0.79, 'Inert'; 'O2', 1, 0.21, 'Oxidant'};
            self.UITable_P.Data = self.UITable_R.Data(:, 1);
        case '3' % AIR (20.9476% O2 + 78.084% N2 + 0.9365% Ar + 0.0319% CO2)
            self.UITable_R.Data = {'N2', 0.78084/0.209476, 0.78084, 'Inert'; 'O2', 1, 0.209476, 'Oxidant'; 'Ar', 0.009365/0.209476, 0.009365, 'Inert'; 'CO2', 0.000319/0.209476, 0.000319, 'Oxidant'};
            self.UITable_P.Data = self.UITable_R.Data(:, 1);
        case '4' % METHANE + AIR IDEAL
            phi_t = 2;
            Ni = [phi_t * 0.79/0.21, phi_t, 1];
            Xi = Ni / sum(Ni);
            self.UITable_R.Data = {'N2', Ni(1), Xi(1), 'Inert'; 'O2', Ni(2), Xi(2), 'Oxidant'; 'CH4', Ni(3), Xi(3), 'Fuel'};
            self.UITable_P.Data = self.UITable_R.Data(:, 1);
        case '5' % METHANE + AIR
            phi_t = 2;
            Ni = [phi_t * 0.78084/0.209476, phi_t, phi_t * 0.009365/0.209476, phi_t * 0.009365/0.209476, 1];
            Xi = Ni / sum(Ni);
            self.UITable_R.Data = {'N2', Ni(1), Xi(1), 'Inert'; 'O2', Ni(2), Xi(2), 'Oxidant'; 'Ar', Ni(3), Xi(3), 'Inert'; 'CO2', Ni(4), Xi(4), 'Oxidant'; 'CH4', Ni(5), Xi(5), 'Fuel'};
            self.UITable_P.Data = self.UITable_R.Data(:, 1);
    end
end