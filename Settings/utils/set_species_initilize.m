function self = set_species_initilize(self, species)
    % Fill the properties matrix with the data of the mixture
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %     species (cell): Species contained in the system
    %
    % Returns:
    %     properties_matrix (float): Properties matrix

    % Get index species
    ind = find_ind(self.S.LS, species);
    % Fill properties matrix
    for i = self.S.NS:-1:1
        self.C.M0.value(ind(i), self.C.M0.ind_W) = self.DB.(species{i}).mm; % [g/mol]
        self.C.M0.value(ind(i), self.C.M0.ind_hfi) = self.DB.(species{i}).hf / 1000; % [kJ/mol]
        self.C.M0.value(ind(i), self.C.M0.ind_efi) = self.DB.(species{i}).ef / 1000; % [kJ/mol]
        self.C.M0.value(ind(i), self.C.M0.ind_phase) = self.DB.(species{i}).phase; % [bool]
    end

end
