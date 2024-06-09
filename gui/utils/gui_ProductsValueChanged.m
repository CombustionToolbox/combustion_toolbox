function gui_ProductsValueChanged(app)
    % Update List of species considered as Products
    try
        if isempty(app.Products.Value)

            if isempty(app.mixture.listSpecies)
                app.listbox_Products.Items = [];
                return
            end
            
            % Get react species
            listSpecies = [app.mixture.listSpeciesFuel, app.mixture.listSpeciesOxidizer];

            % Get list of species that may appear at chemical equilibrium
            app.listbox_Products.Items = findProducts(app.chemicalSystem, listSpecies, 'flag_ion', app.IonizedspeciesCheckBox.Value);
            return
        end

        % Update List of Products depending of the value of the equivalence ratio
        if strcmpi(app.Products.Value, 'Complete Reaction')
            gui_edit_phiValueChanged(app, []); 
            return
        end
        
        FLAG_ADD = ~isempty(combustiontoolbox.utils.findIndex(app.database.listSpecies, app.Products.Value));

        if FLAG_ADD
            app.chemicalSystem.setListSpecies(app.database, app.Products.Value);
            app.listbox_Products.Items = unique([app.listbox_Products.Items, app.chemicalSystem.listSpecies], 'stable');
            return
        end

        app.chemicalSystem.setListSpecies(app.database, app.Products.Value);
        app.listbox_Products.Items = app.chemicalSystem.listSpecies;

    catch
        app.listbox_Products.Items = {};
    end
    
end