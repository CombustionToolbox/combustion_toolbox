function gui_UITable_RCellEdit(obj, event)
    % Update values of the UITable items with the changes made
    try
        % Initialize app (fast: transfer DB)
        app = App('fast', obj.DB_master, obj.DB);
        % Get reactant species
        app = gui_get_reactants(obj, event, app);
        % Compute properties of the mixture
        app = gui_compute_propReactants(obj, app);
        % Update UITable classes
        gui_update_UITable_R(obj, app);
        obj.UITable_P.Data = obj.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
        obj.UITable_R2.Data = obj.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(obj, app);
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end