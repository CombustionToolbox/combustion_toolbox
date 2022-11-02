function properties_matrix = set_species(self, species, moles, T)
    % Fill the properties matrix with the data of the mixture
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %     species (cell): Species contained in the system
    %     moles (float):  Moles of the species in the mixture [mol]
    %     T (float):      Temperature [K]
    %
    % Returns:
    %     properties_matrix (float): Properties matrix

    % Initialize
    properties_matrix = self.C.M0.value;
    % Get index species
    ind = find_ind(self.S.LS, species);
    % Fill properties matrix
    for i = length(moles):-1:1

        if length(self.DB.(species{i}).T) > 1
            h0i = species_h0(species{i}, T, self.DB); % [kJ/mol]
            cPi = species_cP(species{i}, T, self.DB); % [J/mol-K]
            s0i = species_s0(species{i}, T, self.DB); % [kJ/mol-K]
        else
            h0i = self.C.M0.value(ind(i), self.C.M0.ind_hfi); % [kJ/mol]
            cPi = 0; % [J/mol-K]
            s0i = 0; % [kJ/mol-K]
        end

        properties_matrix(ind(i), self.C.M0.ind_ni:end) = [moles(i), h0i, cPi, s0i];
    end

end
