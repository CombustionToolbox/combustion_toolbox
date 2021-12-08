function [obj, temp_app] = gui_edit_phiValueChanged(obj, event)
    % Update moles and mole fractions of the reactant UITables for the
    % given equivalence ratio
    try
        % If UITable_R is empty the equivalence can not change
        if isempty(obj.UITable_R.Data) 
            return
        end
        % Set FLAG compute moles from equivalence ratio
        FLAG_COMPUTE_FROM_PHI = true;
        % Create temporaly app required for preliminary calculations
        temp_app = gui_create_temp_app(obj, event, FLAG_COMPUTE_FROM_PHI);
        % Update UITable classes
        gui_update_UITable_R(obj, temp_app);
        obj.UITable_P.Data = obj.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
        obj.UITable_R2.Data = obj.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F and percentage Fuel
        obj.edit_phi2.Value = obj.edit_phi.Value;
        obj.edit_OF.Value = 1/temp_app.PS.strR{1}.FO;
        obj.edit_F.Value = temp_app.PS.strR{1}.percentage_Fuel;
        % Check List of products (only in case of ListProducts == Complete reaction)
        temp_app = check_ListProducts(obj, temp_app);
    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
    end
end

function app = check_ListProducts(obj, app)
    if strcmpi(obj.Products.Value, 'Complete Reaction')
        if isempty(obj.UITable_R.Data)
            % Set lamp to Error color
            obj.Lamp.Color = obj.color_lamp_warning;
            % Print error
            message = {'First, define initial mixture'};
            uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
        end
        % Update List of Products depending of the value of the equivalence ratio
        phi = app.PS.strR{1}.phi;
        phi_c = app.PS.strR{1}.phi_c;
        app = ListSpecies(obj, obj.Products.Value, phi, phi_c);
        obj.listbox_Products.Items = app.S.LS;
    end
end