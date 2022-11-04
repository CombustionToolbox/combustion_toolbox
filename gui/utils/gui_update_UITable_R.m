function gui_update_UITable_R(app, self)
    % Update data in the UITable_R with the next order: Inert -> Oxidizer -> Fuel
    species = [self.PD.S_Inert, self.PD.S_Oxidizer, self.PD.S_Fuel];
    Nspecies = length(species); 
    if isempty(self.PD.S_Fuel)
        self.PD.N_Fuel = []; % Set to 1 by default
    end
    moles = [self.PD.N_Inert, self.PD.N_Oxidizer, self.PD.N_Fuel];
    molar_fraction = moles/sum(moles); % It is easier to recompute
    typeSpecies = get_typeSpecies(self);
    if ~isempty(app.UITable_R.Data) && Nspecies == length(app.UITable_R.Data(:, 1))
        temperatures = app.UITable_R.Data(:, 5)';
    else
        temperatures = create_cell_ntimes(self.PD.TR.value, Nspecies);
    end
    % Check if is a condensed species with a fixed temperature
    [app, temperatures] = gui_check_temperature_reactants(app, self.DB, species, temperatures, Nspecies);
    % Update table
    app.UITable_R.Data = [species; vector2cell(moles); vector2cell(molar_fraction); typeSpecies; temperatures]';
end