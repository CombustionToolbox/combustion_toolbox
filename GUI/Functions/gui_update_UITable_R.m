function gui_update_UITable_R(obj, app)
    % Update data in the UITable_R with the next order: Inert -> Oxidizer -> Fuel
    species = [app.PD.S_Inert, app.PD.S_Oxidizer, app.PD.S_Fuel];
    Nspecies = length(species); 
    if isempty(app.PD.S_Fuel)
        app.PD.N_Fuel = []; % Set to 1 by default
    end
    moles = [app.PD.N_Inert, app.PD.N_Oxidizer, app.PD.N_Fuel];
    molar_fraction = moles/sum(moles); % It is easier to recompute
    typeSpecies = get_typeSpecies(app);
    if ~isempty(obj.UITable_R.Data) && Nspecies == length(obj.UITable_R.Data(:, 1))
        temperature = obj.UITable_R.Data(:, 5)';
    else
        temperature = create_cell_ntimes(app.PD.TR.value, Nspecies);
    end
    obj.UITable_R.Data = [species; vector2cell(moles); vector2cell(molar_fraction); typeSpecies; temperature]';
end