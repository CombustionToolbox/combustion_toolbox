function gui_update_frozen(app)
    % Set frozen chemistry and update Listbox of species
    if app.FrozenchemistryCheckBox.Value
        try
            species = app.UITable_R.Data(:, 1);
            app.listbox_Products.Items = species;
        catch
            app.listbox_Products.Items = {};
        end
    else
        gui_ProductsValueChanged(app)
    end
end