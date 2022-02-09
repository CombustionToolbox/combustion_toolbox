function gui_ProductsValueChanged(obj)
    % Update List of species considered as Products
    try
        if isempty(obj.Products.Value)
            temp_app = ListSpecies(obj);
            obj.listbox_Products.Items = temp_app.S.LS;
        else
            % Update List of Products depending of the value of the equivalence ratio
            if strcmpi(obj.Products.Value, 'Complete Reaction')
                gui_edit_phiValueChanged(obj, []); 
            else
                temp_app = ListSpecies(obj, obj.Products.Value);
                obj.listbox_Products.Items = [obj.listbox_Products.Items, temp_app.S.LS];
            end
        end
    catch
        obj.listbox_Products.Items = {};
    end
end