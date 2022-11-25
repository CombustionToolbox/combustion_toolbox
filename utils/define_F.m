function self = define_F(self)
    % Set Fuel of the mixture
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    if ~isempty(self.PD.S_Fuel)
        self.PD.R_Fuel = set_species(self, self.PD.S_Fuel, self.PD.N_Fuel, self.PD.TR.value);
        self.PS.strR_Fuel = compute_properties(self, self.PD.R_Fuel, self.PD.pR.value, self.PD.TR.value);
        self.PS.strR_Fuel = assign_values_elements(self, self.PS.strR_Fuel);
        self.PD.Fuel.x = self.PS.strR_Fuel.x;
        self.PD.Fuel.x2 = self.PS.strR_Fuel.x2;
        self.PD.Fuel.x3 = self.PS.strR_Fuel.x3;
        self.PD.Fuel.y = self.PS.strR_Fuel.y;
        self.PD.Fuel.z = self.PS.strR_Fuel.z;
        self.PD.Fuel.w = self.PS.strR_Fuel.w;
        self.PD.phi_t = abs(self.PD.Fuel.x + self.PD.Fuel.x2 + ...
            self.PD.Fuel.x3 + self.PD.Fuel.y / 4 + ...
            - self.PD.Fuel.z / 2) / (0.5 * self.S.atoms_ox_ref);
    else
        self.PD.R_Fuel = 0; self.PD.phi_t = 1;
        self.PD.Fuel.x = 0;
        self.PD.Fuel.x2 = 0;
        self.PD.Fuel.x3 = 0;
        self.PD.Fuel.y = 0;
        self.PD.Fuel.z = 0;
        self.PD.Fuel.w = 0;
        self.C.FLAG_Fuel = 0;
    end

end

% SUB-PASS FUNCTIONS
function mix = assign_values_elements(self, mix)
    % Assign values for C, H, O, and N elements

    if isempty(self.E.ind_C), mix.x = 0; else, mix.x = mix.NatomE(self.E.ind_C); end
    if isempty(self.E.ind_H), mix.y = 0; else, mix.y = mix.NatomE(self.E.ind_H); end
    if isempty(self.E.ind_O), mix.z = 0; else, mix.z = mix.NatomE(self.E.ind_O); end
    if isempty(self.E.ind_N), mix.w = 0; else, mix.w = mix.NatomE(self.E.ind_N); end
    if isempty(self.E.ind_S), mix.x2 = 0; else, mix.x2 = mix.NatomE(self.E.ind_S); end
    if isempty(self.E.ind_Si), mix.x3 = 0; else, mix.x3 = mix.NatomE(self.E.ind_Si); end
end
