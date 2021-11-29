function gui_UITable_RCellEdit(obj, event)
    % Update values of the UITable items with the changes made
    try
        % Initialize app (fast: transfer DB)
        app = App('fast', obj.DB_master, obj.DB);
        % Get reactant species
        app = gui_get_reactants(obj, event, app);
        % Compute properties of the mixture
        app = gui_compute_propReactants(obj, app);
        % Update temperature of the mixture (if needed)
        update_temperature_mixture(obj, app);
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

function update_temperature_mixture(obj, app)
    % Get temperature of the mixture (if needed)
    temperature = compute_temperature_mixture(app, obj.UITable_R.Data(:, 1), obj.UITable_R.Data(:, 2), obj.UITable_R.Data(:, 5));
    obj.PR1.Value = sprintf('%.4g', temperature);
end