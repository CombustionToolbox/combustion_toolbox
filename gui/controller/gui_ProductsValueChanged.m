function gui_ProductsValueChanged(obj)
    % Update List of species considered as Products
    try
        if isempty(obj.Products.Value)
            temp_app = list_species(obj);
            obj.listbox_Products.Items = temp_app.S.LS;
        else
            % Update List of Products depending of the value of the equivalence ratio
            if strcmpi(obj.Products.Value, 'Complete Reaction')
                gui_edit_phiValueChanged(obj, []); 
            else
                try
                    temp_app = list_species(obj, obj.Products.Value);
                    temp_app.DB.(obj.Products.Value);
                    obj.listbox_Products.Items = unique([obj.listbox_Products.Items, temp_app.S.LS], 'stable');
                catch
                    temp_app = list_species(obj, obj.Products.Value);
                    obj.listbox_Products.Items = temp_app.S.LS;
                end
            end
        end
    catch
        obj.listbox_Products.Items = {};
    end
end