function gui_UITable_RCellEdit(app, event)
    % Update values of the UITable items with the changes made
    
    try
        % Definitions
        temperature = gui_get_prop(app.PR1.Value, 'first'); % [K]
        pressure = gui_get_prop(app.PR2.Value, 'first');    % [bar]
        listSpecies = app.listbox_Products.Items;
        
        % Define chemical system
        app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database, listSpecies);

        % Get previous ratioOxidizer
        ratioOxidizer = app.mixture.ratioOxidizer;
        
        % Initialize mixture
        app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);

        % Define chemical state
        gui_get_reactants(app, event);
        
        % Update temperature of the mixture (if needed)
        temperature = updateTemperature(app);

        % Define properties
        app.mixture = setProperties(app.mixture, 'temperature', temperature, 'pressure', pressure);
        
        % Set ratioOxidizer
        app.mixture.ratioOxidizer = ratioOxidizer;
        
        % Update UITable classes
        gui_update_UITable_R(app);
        app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)

        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(app);

        % Update GUI: Listbox species considered as products.
        % In case frozen chemistry or ionization is considered.
        gui_update_frozen(app);




        % % Initialize self (fast: transfer DB)
        % self = App('fast', app.DB_master, app.DB, app.S.LS);
        % % Get reactant species
        % self = gui_get_reactants(app, event, self);
        % % Compute properties of the mixture
        % self = gui_compute_propReactants(app, self);
        % % Update temperature of the mixture (if needed)
        % update_temperature_mixture(app, self);
        % % Update UITable classes
        % gui_update_UITable_R(app, self);
        % app.UITable_P.Data = app.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
        % app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % % Update GUI: equivalence ratio, O/F and percentage Fuel
        % gui_update_phi(app, self);
    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
    end

end

function temperature = updateTemperature(app)
    % Get temperature of the mixture (if needed)
    temperature = compute_temperature_mixture(app.ProblemType.Value, app.chemicalSystem.species, app.UITable_R.Data(:, 1), app.UITable_R.Data(:, 2), app.UITable_R.Data(:, 5));
    app.PR1.Value = sprintf('%.4g', temperature);
end