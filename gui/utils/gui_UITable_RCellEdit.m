function gui_UITable_RCellEdit(app, event)
    % Update values of the UITable items with the changes made
    
    try
        % Definitions
        pressure = gui_get_prop(app.PR2.Value, 'first'); % [bar]
        listSpecies = app.listbox_Products.Items;
        
        % Define chemical system
        if isempty(listSpecies)
            app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database);
        else
            app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database, listSpecies);
        end
        
        % Initialize mixture
        app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);

        % Define chemical state
        gui_get_reactants(app, event);
        
        % Set temperature of the mixture
        setTemperature(app);

        % Set pressure of the mixture
        setPressure(app.mixture, pressure);
        
        % Set ratioOxidizer based on current oxidizer species
        if ~isempty(app.mixture.stoichiometricMoles)
            app.mixture.ratioOxidizer = app.mixture.ratioOxidizer / app.mixture.stoichiometricMoles * app.mixture.equivalenceRatio;
        end
        
        % Update UITable classes
        gui_update_UITable_R(app);
        app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)

        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(app);

        % Update GUI: Listbox species considered as products.
        % In case frozen chemistry or ionization is considered.
        gui_update_frozen(app);

    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
    end

end

function setTemperature(app)
    % Set temperature of the mixture

    % Definitions
    app.mixture.problemType = app.ProblemType.Value;
    speciesTemperatures = cell2vector(app.UITable_R.Data(:, 5))';

    % Reorganize species temperature in the same order as in listSpecies
    [~, index] = ismember(app.mixture.listSpecies, app.UITable_R.Data(:, 1));
    speciesTemperatures = speciesTemperatures(index);

    % Compute equilibrium temperature
    setTemperatureSpecies(app.mixture, speciesTemperatures);

    % Update GUI
    app.PR1.Value = sprintf('%.4g', app.mixture.T);
end