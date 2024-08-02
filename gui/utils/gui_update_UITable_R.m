function gui_update_UITable_R(app)
    % Update data in the UITable_R with the next order: Inert -> Oxidizer -> Fuel

    % Get chemical species and molar composition
    listSpecies = [app.mixture.listSpeciesInert, app.mixture.listSpeciesOxidizer, app.mixture.listSpeciesFuel];
    numSpecies = length(listSpecies);
    moles = [app.mixture.molesInert, app.mixture.molesOxidizer, app.mixture.molesFuel];
    molar_fraction = moles / sum(moles);

    % Catalog chemical species (fuel, oxidizer, or inert)
    typeSpecies = getTypeSpecies(app.mixture);

    if ~isempty(app.UITable_R.Data) && numSpecies == length(app.UITable_R.Data(:, 1))
        temperatures = app.UITable_R.Data(:, 5)';
    else
        temperatures = repmat({app.mixture.T}, [1, numSpecies]);
    end

    % Check if is a condensed species with a fixed temperatures
    [app, temperatures] = gui_check_temperature_reactants(app, listSpecies, temperatures, numSpecies);

    % Update table
    app.UITable_R.Data = [listSpecies; vector2cell(moles); vector2cell(molar_fraction); typeSpecies; temperatures]';
end