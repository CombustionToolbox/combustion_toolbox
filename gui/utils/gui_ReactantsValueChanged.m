function gui_ReactantsValueChanged(app, event)
    % Update values of the UITable items with:
    % - a given predefined set of reactants
    % - the new species added in the finder
    try
        % Initialize self (fast: transfer DB)
        self = App('fast', app.DB_master, app.DB);
        % If the empty item was selected clear UITable and return
        if app.Reactants.Value == 1 % No species selected
            gui_empty_UITables(app);
            return
        elseif isempty(app.Reactants.Value)
            return
        end
        % Set reactant species
        self = gui_set_reactants(app, event, self);
        % Compute properties of the mixture
        self = gui_compute_propReactants(app, self);
        % Update UITable classes
        gui_update_UITable_R(app, self);
        app.UITable_R2.Data = app.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
        % Update GUI: equivalence ratio, O/F and percentage Fuel
        gui_update_phi(app, self);
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
function self = gui_set_reactants(obj, event, self)
    % Set reactant species from the standard set drop down

    % Check type of air
    % - ideal:     21.0000% O2, 79.0000% N2
    % - non-ideal: 21.0000% O2, 78.0840% N2, 0.9365% Ar, 0.0319% CO2
    FLAG_IDEAL_AIR = obj.IdealAirCheckBox.Value;
    % Get equivalence ratio (phi) from GUI
    if ~strcmp(obj.edit_phi.Value, '-') % Default value phi = 1 (stoichiometric)
        self.PD.phi.value = gui_get_prop(self, 'phi', obj.edit_phi.Value);
    end
    % Set species in the mixture
    switch obj.Reactants.Value
        case 2 % Air
            self = set_air(self, FLAG_IDEAL_AIR);
        case 3 % Methane + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'CH4'};
        case 4 % Ethane + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C2H6'};
        case 5 % Propane + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C3H8'};
        case 6 % Acetylene + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C2H2_acetylene'};
        case 7 % Ethylene + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C2H4'};
        case 8 % Bencene + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C6H6'};
        case 9 % Iso-octane + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C8H18_isooctane'};
        case 10 % Hydrogen + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'H2'};
        case 11 % Carbon monoxide + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'CO'};
        case 12 % Methanol + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'CH3OH'};
        case 13 % Ethanol + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'C2H5OH'};
        case 14 % Natural Gas + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'CH4','C2H6','C3H8'};
            self.PD.N_Fuel = [0.85, 0.1, 0.05];
        case 15 % Syngas + Air
            self = set_air(self, FLAG_IDEAL_AIR);
            self.PD.S_Fuel = {'CO','H2'};  
            self.PD.N_Fuel = [0.5, 0.5];
        case 16 % LH2 + LOX
            self.PD.S_Fuel = {'H2bLb'};
            self.PD.S_Oxidizer = {'O2bLb'};
            obj.Products.Value = 'Hydrogen (L)';
            % Update List of species considered as Products
            gui_ProductsValueChanged(obj);
            % Update Listbox (extended settings)
            obj.listbox_LS.Items = obj.listbox_Products.Items;
        case 17 % RP1 + LOX
            self.PD.S_Fuel = {'RP_1'};
            self.PD.S_Oxidizer = {'O2bLb'};
        otherwise % SET NEW SPECIES
            try
                species = gui_seeker_exact_value(obj, event, self.S.LS_DB);
            catch
                message = {'Species not found.'};
                uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
                return
            end
            % Get data of the current mixture
            if ~isempty(obj.UITable_R.Data)
                % Check if the species is already included          
                self = gui_get_reactants(obj, event, self);
            end
            % Add new species to the mixture (fuel by default)
            if ~any(find_ind(self.PD.S_Fuel, species))
                self.PD.S_Fuel = [self.PD.S_Fuel, {species}];
                self.PD.N_Fuel = [self.PD.N_Fuel, 1];
            end
    end
end