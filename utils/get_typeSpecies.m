function typeSpecies = get_typeSpecies(mix)
    % Create cell array with the type of species in the mixture
    %
    % Args:
    %     mix (Mixture): Mixture class
    %
    % Returns:
    %     typeSpecies (cell): Cell array with the type of species in the mixture

    typeFuel = create_cell_ntimes('Fuel', length(mix.listSpeciesFuel));
    typeOxidizer = create_cell_ntimes('Oxidizer', length(mix.listSpeciesOxidizer));
    typeInert = create_cell_ntimes('Inert', length(mix.listSpeciesInert));
    typeSpecies = [typeInert, typeOxidizer, typeFuel];
end