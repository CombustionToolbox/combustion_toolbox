function gui_update_frozen(obj)
    % Set frozen chemistry and update Listbox of species
    if obj.FrozenchemistryCheckBox.Value
        try
            species = obj.UITable_R.Data(:, 1);
            obj.listbox_Products.Items = species;
        catch
            obj.listbox_Products.Items = {};
        end
    else
        gui_ReactionValueChanged(obj)
    end
end