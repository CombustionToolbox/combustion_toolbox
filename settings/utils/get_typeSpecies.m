function typeSpecies = get_typeSpecies(self)
    % Create cell array with the type of species in the mixture
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    typeFuel = create_cell_ntimes('Fuel', length(self.PD.N_Fuel));
    typeOxidizer = create_cell_ntimes('Oxidizer', length(self.PD.N_Oxidizer));
    typeInert = create_cell_ntimes('Inert', length(self.PD.N_Inert));
    typeSpecies = [typeInert, typeOxidizer, typeFuel];
end