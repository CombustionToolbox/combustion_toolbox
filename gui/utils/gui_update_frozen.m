function gui_update_frozen(app)
    % Set frozen chemistry and update Listbox of species
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    
    % Check if frozen chemistry is selected
    if ~app.FrozenchemistryCheckBox.Value
        gui_ProductsValueChanged(app)

        % Update Listbox (extended settings)
        public_ProductsValueChanged(app);
        return
    end
    
    if contains(app.ProblemType.Value, 'ROCKET', 'IgnoreCase', true)
        return
    end
    
    % Set list species as list of reactants
    try
        species = app.UITable_R.Data(:, 1);
        app.listbox_Products.Items = species;
    catch
        app.listbox_Products.Items = {};
    end

    % Update Listbox (extended settings)
    public_ProductsValueChanged(app);
end