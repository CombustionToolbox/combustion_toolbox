function gui_ProductsValueChanged(app)
    % Update List of species considered as Products
    try
        if isempty(app.Products.Value)
            temp_app = list_species(app);
            app.listbox_Products.Items = temp_app.S.LS;
        else
            % Update List of Products depending of the value of the equivalence ratio
            if strcmpi(app.Products.Value, 'Complete Reaction')
                gui_edit_phiValueChanged(app, []); 
            else
                try
                    temp_app = list_species(app, app.Products.Value);
                    temp_app.DB.(app.Products.Value);
                    app.listbox_Products.Items = unique([app.listbox_Products.Items, temp_app.S.LS], 'stable');
                catch
                    temp_app = list_species(app, app.Products.Value);
                    app.listbox_Products.Items = temp_app.S.LS;
                end
            end
        end
    catch
        app.listbox_Products.Items = {};
    end
end