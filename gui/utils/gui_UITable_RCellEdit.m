function gui_UITable_RCellEdit(app, event)
    % Update values of the UITable items with the changes made
    
    try
        % Initialize self (fast: transfer DB)
        self = App('fast', app.DB_master, app.DB);
        % Get reactant species
        self = gui_get_reactants(app, event, self);
        % Compute properties of the mixture
        self = gui_compute_propReactants(app, self);
        % Update temperature of the mixture (if needed)
        update_temperature_mixture(app, self);
        % Update UITable classes
        gui_update_UITable_R(app, self);
        app.UITable_P.Data = app.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
        app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(app, self);
    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
    end
end

function update_temperature_mixture(obj, self)
    % Get temperature of the mixture (if needed)
    temperature = compute_temperature_mixture(self, obj.UITable_R.Data(:, 1), obj.UITable_R.Data(:, 2), obj.UITable_R.Data(:, 5));
    obj.PR1.Value = sprintf('%.4g', temperature);
end