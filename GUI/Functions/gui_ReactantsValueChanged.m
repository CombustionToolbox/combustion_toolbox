function gui_ReactantsValueChanged(obj, event)
    % Update values of the UITable items with:
    % - a given predefined set of reactants
    % - the new species added in the finder
    try
        % Initialize app (fast: transfer DB)
        app = App('fast', obj.DB_master, obj.DB);
        % If the empty item was selected clear UITable and return
        if strcmp(obj.Reactants.Value, '1') % No species selected
            gui_empty_UITables(obj);
            return
        elseif isempty(obj.Reactants.Value)
            return
        end
        % Set reactant species
        app = gui_set_reactants(obj, event, app);
        % Compute properties of the mixture
        app = gui_compute_propReactants(obj, app);
        % Update UITable classes
        gui_update_UITable_R(obj, app);
        obj.UITable_R2.Data = obj.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(obj, app);
        % Update GUI: Listbox species considered as products.
        % In case frozen chemistry or ionization is considered.
        gui_update_frozen(obj);
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end

% SUB-PASS FUNCTIONS
function app = gui_set_reactants(obj, event, app)
    % Set reactant species from the standard set drop down

    % Check type of air
    % - ideal:     21.0000% O2, 79.0000% N2
    % - non-ideal: 21.0000% O2, 78.0840% N2, 0.9365% Ar, 0.0319% CO2
    FLAG_IDEAL_AIR = obj.IdealAirCheckBox.Value;
    % Get equivalence ratio (phi) from GUI
    if ~strcmp(obj.edit_phi.Value, '-') % Default value phi = 1 (stoichiometric)
        app.PD.phi.value = gui_get_prop(app, 'phi', obj.edit_phi.Value);
    end
    % Set species in the mixture
    switch obj.Reactants.Value
        case '2' % Air
            app = set_air(app, FLAG_IDEAL_AIR);
        case '3' % Methane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CH4'};
        case '4' % Ethane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H6'};
        case '5' % Propane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C3H8'};
        case '6' % Acetylene + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H2_acetylene'};
        case '7' % Ethylene + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H4'};
        case '8' % Bencene + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C6H6'};
        case '9' % Iso-octane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C8H18_isooctane'};
        case '10' % Hydrogen + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'H2'};
        case '11' % Carbon monoxide + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CO'};
        case '12' % Methanol + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CH3OH'};
        case '13' % Ethanol + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H5OH'};
        case '14' % Natural Gas + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CH4','C2H6','C3H8'};
            app.PD.N_Fuel = [0.85, 0.1, 0.05];
        case '15' % Syngas + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CO','H2'};  
            app.PD.N_Fuel = [0.5, 0.5];
        otherwise % SET NEW SPECIES
            try
                species = gui_seeker_species(obj, event);
            catch
                message = {'Species not found.'};
                uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
                return
            end
            % Get data of the current mixture
            if ~isempty(obj.UITable_R.Data)
                app = gui_get_reactants(obj, event, app);
            end
            % Add new species to the mixture (fuel by default)
            app.PD.S_Fuel = [app.PD.S_Fuel, {species}];
            app.PD.N_Fuel = [app.PD.N_Fuel, 0];
    end
end