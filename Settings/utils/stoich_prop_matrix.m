function self = stoich_prop_matrix(self)
    % Initialize the stoichiometric matrix and properties matrix
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct):  Data of the mixture, conditions, and databases

    self.C.A0.value = zeros(self.S.NS, self.E.NE);
    self.C.M0.value = zeros(self.S.NS, self.C.N_prop.value);

    for i = 1:self.S.NS
        txFormula = self.DB.(self.S.LS{i}).txFormula;
        self.DB.(self.S.LS{i}).Element_matrix = set_element_matrix(txFormula, self.E.elements);
        self.C.A0.value(i, self.DB.(self.S.LS{i}).Element_matrix(1, :)) = self.DB.(self.S.LS{i}).Element_matrix(2, :);
    end

    % Fill constant values of properties matrix
    self = set_species_initilize(self, self.S.LS);
    % Initialize vector [moles, phase]
    self.C.N0.value = self.C.M0.value(:, [self.C.M0.ind_ni, self.C.M0.ind_phase]);
end
