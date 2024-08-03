function gui_update_ions(app)
    % Set ionized species and update Listbox of species
    %
    % Args:
    %     app (object): Combustion Toolbox app object

    % Set product list to empty (auto-selection)
    app.Products.Value = [];
    
    % Get react species
    listSpecies = [app.mixture.listSpeciesFuel, app.mixture.listSpeciesOxidizer];

    % Get list of species that may appear at chemical equilibrium
    app.listbox_Products.Items = findProducts(app.chemicalSystem, listSpecies, 'flag_ion', app.IonizedspeciesCheckBox.Value);

    % Update Listbox (extended settings)
    public_ProductsValueChanged(app);
end