function gui_ReactantsValueChanged(app, event)
    % Update values of the UITable items with:
    %   * a given predefined set of reactants
    %   * the new species added in the finder
    try
        % If the empty item was selected clear UITable and return
        if app.Reactants.Value == 1 % No species selected
            gui_empty_UITables(app);
            return
        elseif isempty(app.Reactants.Value)
            return
        end
        
        % Definitions
        temperature = gui_get_prop(app.PR1.Value, 'first'); % [K]
        pressure = gui_get_prop(app.PR2.Value, 'first');    % [bar]
        listSpecies = app.listbox_Products.Items;

        % Define chemical system
        if isempty(listSpecies)
            app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database);
        else
            app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database, listSpecies);
        end
        
        % Temp
        ratioOxidizer = app.mixture.ratioOxidizer;
        stoichiometricMoles = app.mixture.stoichiometricMoles;

        % Initialize mixture
        app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);
        if ~isempty(stoichiometricMoles)
            app.mixture.ratioOxidizer = ratioOxidizer;
        end

        % Set chemical state
        gui_set_reactants(app, event);

        % Define properties
        if strcmp(app.edit_phi.Value, '-')
            app.mixture = setProperties(app.mixture, 'temperature', temperature, 'pressure', pressure);
        else
            equivalenceRatio = gui_get_prop(app.edit_phi.Value, 'first'); % [-]
            app.mixture = setProperties(app.mixture, 'temperature', temperature, 'pressure', pressure, 'equivalenceRatio', equivalenceRatio);
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
        message = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
        fprintf('%s\n', message);
    end
    
end

% SUB-PASS FUNCTIONS
function app = gui_set_reactants(app, event)
    % Set reactant species from the standard set drop down
    
    % Import packages
    import combustiontoolbox.utils.getAir
    import combustiontoolbox.utils.findIndex

    % Definitions
    FLAG_FUEL = true;
    FLAG_OXIDIZER = false;

    % Check type of air
    % - ideal:     21.0000% O2, 79.0000% N2
    % - non-ideal: 21.0000% O2, 78.0840% N2, 0.9365% Ar, 0.0319% CO2
    FLAG_IDEAL_AIR = app.IdealAirCheckBox.Value;
    
    % Set species in the mixture
    switch app.Reactants.Value
        case 2 % Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR);
            FLAG_FUEL = false; FLAG_OXIDIZER = true; app.edit_phi.Value = '-';
        case 3 % Methane + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'CH4'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 4 % Ethane + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C2H6'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 5 % Propane + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C3H8'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 6 % Acetylene + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C2H2_acetylene'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 7 % Ethylene + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C2H4'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 8 % Bencene + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C6H6'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 9 % Iso-octane + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C8H18_isooctane'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 10 % Hydrogen + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'H2'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 11 % Carbon monoxide + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'CO'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 12 % Methanol + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'CH3OH'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 13 % Ethanol + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'C2H5OH'}; molesFuel = 1; app.edit_phi.Value = '1';
        case 14 % Natural Gas + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'CH4','C2H6','C3H8'}; molesFuel = [0.85, 0.1, 0.05]; app.edit_phi.Value = '1';
        case 15 % Syngas + Air
            [listSpeciesOxidizer, molesOxidizer] = getAir(FLAG_IDEAL_AIR); FLAG_OXIDIZER = true;
            listSpeciesFuel = {'CO','H2'}; molesFuel = [0.5, 0.5]; app.edit_phi.Value = '1';
        case 16 % LH2 + LOX
            listSpeciesFuel = {'H2bLb'}; molesFuel = 1; app.edit_phi.Value = '1';
            listSpeciesOxidizer = {'O2bLb'}; molesOxidizer = 1; FLAG_OXIDIZER = true;
        case 17 % RP1 + LOX
            listSpeciesFuel = {'RP_1'}; molesFuel = 1; app.edit_phi.Value = '1';
            listSpeciesOxidizer = {'O2bLb'}; molesOxidizer = 1; FLAG_OXIDIZER = true;
        otherwise % SET NEW SPECIES

            try
                listSpeciesAdd = gui_seeker_exact_value(app, event, app.database.listSpecies);
            catch
                message = {'Species not found.'};
                uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
                return
            end

            % Get data of the current mixture
            if ~isempty(app.UITable_R.Data)
                gui_get_reactants(app, event);
            end

            % Add new species to the mixture (inert by default)
            if any(findIndex(app.mixture.listSpecies, listSpeciesAdd))
                return
            end

            molesAdd = 1;
            set(app.mixture, {listSpeciesAdd}, 'inert', molesAdd);
            % app.mixture.ratioOxidizer = molesOxidizer;
            return
    end
    
    % Set chemical state
    if FLAG_FUEL
        set(app.mixture, listSpeciesFuel, 'fuel', molesFuel);
    end

    if FLAG_OXIDIZER
        set(app.mixture, listSpeciesOxidizer, 'oxidizer', molesOxidizer);
        % app.mixture.ratioOxidizer = molesOxidizer;
    end
    
end