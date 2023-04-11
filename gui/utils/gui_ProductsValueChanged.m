function gui_ProductsValueChanged(app)
    % Update List of species considered as Products
    try
        if isempty(app.Products.Value)
            temp_app = list_species(app);
            app.listbox_Products.Items = temp_app.S.LS;

            % This will be included in the new release
            % species = app.LS_reactants;
            % temp_app = App('fast', app.DB_master, app.DB);
            % LS = find_products(temp_app, species);
            % app.listbox_Products.Items = LS;
            return
        end

        % Update List of Products depending of the value of the equivalence ratio
        if strcmpi(app.Products.Value, 'Complete Reaction')
            gui_edit_phiValueChanged(app, []); 
            return
        end
    
        try
            temp_app = list_species(app, app.Products.Value);
            temp_app.DB.(app.Products.Value);
            app.listbox_Products.Items = unique([app.listbox_Products.Items, temp_app.S.LS], 'stable');
        catch
            temp_app = list_species(app, app.Products.Value);
            app.listbox_Products.Items = temp_app.S.LS;
        end

    catch
        app.listbox_Products.Items = {};
    end
    
end