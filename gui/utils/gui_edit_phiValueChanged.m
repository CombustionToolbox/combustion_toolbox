function [app, temp_app] = gui_edit_phiValueChanged(app, event)
    % Update moles and mole fractions of the reactant UITables for the
    % given equivalence ratio
    try
        % If UITable_R is empty the equivalence can not change
        if isempty(app.UITable_R.Data) 
            return
        end
        % Set FLAG compute moles from equivalence ratio
        FLAG_COMPUTE_FROM_PHI = true;
        % Create temporaly app required for preliminary calculations
        temp_app = gui_create_temp_app(app, event, FLAG_COMPUTE_FROM_PHI);
        % Update UITable classes
        gui_update_UITable_R(app, temp_app);
        app.UITable_P.Data = app.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
        app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F, percentage Fuel, and
        % temperature of reactants
        app.edit_phi2.Value = app.edit_phi.Value;
        app.edit_phi3.Value = app.edit_phi.Value;
        app.edit_OF.Value = temp_app.PS.strR{1}.OF;
        app.edit_F.Value = temp_app.PS.strR{1}.percentage_Fuel;
        app.PR1.Value = sprintf('%g', round(temp_app.PD.TR.value, 2));
        % Check List of products (only in case of ListProducts == Complete reaction)
        temp_app = check_ListProducts(app, temp_app);
    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
    end
end

function self = check_ListProducts(app, self)
    if strcmpi(app.Products.Value, 'Complete Reaction')
        if isempty(app.UITable_R.Data)
            % Set lamp to Error color
            app.Lamp.Color = app.color_lamp_warning;
            % Print error
            message = {'First, define initial mixture'};
            uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
        end
        % Update List of Products depending of the value of the equivalence ratio
        self = list_species(app, app.Products.Value);
        app.listbox_Products.Items = self.S.LS;
    end
end