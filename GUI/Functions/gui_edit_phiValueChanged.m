function gui_edit_phiValueChanged(obj, event)
    % Update moles and mole fractions of the reactant UITables for the
    % given equivalence ratio
    try
        % If UITable_R is empty the equivalence can not change
        if isempty(obj.UITable_R.Data) 
            return
        end
        % Initialize app (fast: transfer DB)
        app = App('fast', obj.DB_master, obj.DB);
        % Set FLAG compute moles from equivalence ratio
        FLAG_COMPUTE_FROM_PHI = true;
        % Get reactant species
        app = gui_get_reactants(obj, event, app, FLAG_COMPUTE_FROM_PHI);
        % Compute properties of the mixture
        app = gui_compute_propReactants(obj, app);
        % Update UITable classes
        gui_update_UITable_R(obj, app);
        obj.UITable_P.Data = obj.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
        obj.UITable_R2.Data = obj.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(obj, app);
    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
    end
end