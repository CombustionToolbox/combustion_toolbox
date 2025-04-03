function gui_edit_phiValueChanged(app, event)
    % Update moles and mole fractions of the reactant UITables for the
    % given equivalence ratio
    try
        if strcmp(app.edit_phi.Value, '-')
            return
        end

        % If UITable_R is empty the equivalence can not change
        if isempty(app.UITable_R.Data)
            % Set lamp to Error color
            app.Lamp.Color = app.color_lamp_warning;
            % Print error
            message = {'First, define initial mixture'};
            uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
            return
        end
        
        % Definitions
        pressure = gui_get_prop(app.PR2.Value, 'first');    % [bar]
        equivalenceRatio = gui_get_prop(app.edit_phi.Value, 'first'); % [-]
        
        % Temporal mixture
        tempMixture = app.mixture.copy();
    
        % Get list species
        if strcmpi(app.Products.Value, 'complete reaction')
            listSpecies = 'complete';
        else
            listSpecies = app.UITable_R.Data(:, 1)';
        end

        % Define chemical system
        app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database, listSpecies);

        % Initialize mixture
        app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);
        
        % Define chemical state
        set(app.mixture, tempMixture.listSpeciesFuel, 'fuel', tempMixture.molesFuel);
        set(app.mixture, tempMixture.listSpeciesOxidizer, 'oxidizer', tempMixture.ratioOxidizer);
        set(app.mixture, tempMixture.listSpeciesInert, 'inert', tempMixture.molesInert);

        % Set temperature of the mixture
        setTemperature(app);

        % Define properties
        app.mixture = setProperties(app.mixture, 'pressure', pressure, 'equivalenceRatio', equivalenceRatio);
        
        % Update UITable classes
        gui_update_UITable_R(app);
        app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, number of moles, mole fractions)

        % Update GUI: equivalence ratio, O/F and percentage Fuel
        app.edit_phi2.Value = app.edit_phi.Value;
        app.edit_phi3.Value = app.edit_phi.Value;
        app.edit_OF.Value = app.mixture.oxidizerFuelMassRatio;
        app.edit_F.Value = app.mixture.percentageFuel;

        % Update GUI: Listbox species considered as products.
        % In case frozen chemistry or ionization is considered.
        % gui_update_frozen(app);

        % Check List of products (only in case of ListProducts == Complete reaction)
        check_ListProducts(app, equivalenceRatio);

    catch ME
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
    end

end

% SUB-PASS FUNCTIONS
function check_ListProducts(app, equivalenceRatio)
    % Check list of products if we have a complete reaction
    if ~strcmpi(app.Products.Value, 'Complete Reaction')
        return
    end

    % Update List of Products depending of the value of the equivalence ratio
    app.chemicalSystem.setListSpecies('complete', equivalenceRatio, app.mixture.equivalenceRatioSoot);
    app.listbox_Products.Items = app.chemicalSystem.listSpecies;
end

function setTemperature(app)
    % Set temperature of the mixture
    
    % Definitions
    app.mixture.problemType = app.ProblemType.Value;
    speciesTemperatures = cell2vector(app.UITable_R.Data(:, 5))';

    % Update table compostion
    gui_update_UITable_R(app);

    % Reorganize species temperature in the same order as in listSpecies
    [~, index] = ismember(app.mixture.listSpecies, app.UITable_R.Data(:, 1));
    speciesTemperatures = speciesTemperatures(index);

    % Compute equilibrium temperature
    setTemperatureSpecies(app.mixture, speciesTemperatures);
    
    % Update GUI
    app.PR1.Value = sprintf('%.4g', app.mixture.T);
end