function typeSpecies = getTypeSpecies(obj)
    % Create cell array with the type of species in the mixture
    %
    % Args:
    %     obj (Mixture): Mixture class
    %
    % Returns:
    %     typeSpecies (cell): Cell array with the type of species in the mixture

    typeFuel = repmat({'Fuel'}, size(obj.listSpeciesFuel));
    typeOxidizer = repmat({'Oxidizer'}, size(obj.listSpeciesOxidizer));
    typeInert = repmat({'Inert'}, size(obj.listSpeciesInert));
    typeSpecies = [typeInert, typeOxidizer, typeFuel];
end