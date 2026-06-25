classdef AppSpeciesPanel < handle
    % Manages GUI species setup, mixture composition, products, and equivalence ratio
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %
    % Examples:
    %     * panel = combustiontoolbox.gui.AppSpeciesPanel(app);
    %     * panel.onReactantsChanged(event);

    properties (Access = private)
        app
        session
    end

    methods
        function obj = AppSpeciesPanel(app, session)
            % AppSpeciesPanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %     session (AppSession): Long-lived GUI services
            %
            % Returns:
            %     obj (AppSpeciesPanel): Initialized species panel object
            if nargin < 1
                app = [];
            end

            if nargin < 2
                session = [];
            end

            obj.app = app;
            obj.session = session;
        end

        function onReactantsChanged(obj, event)
            % Handle predefined or custom reactant selection
            %
            % Args:
            %     event (object): App Designer callback event
            value = obj.componentValue('Reactants', []);

            if isempty(value)
                return
            end

            if isnumeric(value) && value == 1
                obj.clearReactantsTables();
                return
            end

            obj.rebuildChemicalSystem(obj.productSpecies());
            previousRatioOxidizer = obj.app.mixture.ratioOxidizer;
            previousStoichiometricMoles = obj.app.mixture.stoichiometricMoles;
            obj.app.mixture = combustiontoolbox.core.Mixture(obj.app.chemicalSystem);

            if ~isempty(previousStoichiometricMoles)
                obj.app.mixture.ratioOxidizer = previousRatioOxidizer;
            end

            obj.applyReactantsSelection(value, event);
            obj.applyMixturePropertiesFromInputs();
            obj.updateReactantsTable();
            obj.updatePhiControls();
            obj.updateFrozenProducts();
        end

        function onProductsChanged(obj)
            % Handle product preset or custom species selection
            obj.updateProductSpecies();
            obj.updateProductDisplayLists();
        end

        function setCustomProductSpecies(obj, species)
            % Apply product species selected outside the products dropdown
            if isempty(species)
                species = {};
            end

            obj.setProductItems(unique(species, 'stable'));
            obj.updateProductDisplayLists();
        end

        function onReactantsTableEdited(obj, ~)
            % Rebuild mixture from edited reactants table
            if isempty(obj.componentData('UITable_R', {}))
                return
            end

            pressure = obj.parseComponentValue('PR2', 'first');
            obj.rebuildChemicalSystem(obj.productSpecies());
            obj.app.mixture = combustiontoolbox.core.Mixture(obj.app.chemicalSystem);
            obj.setMixtureFromReactantsTable();
            obj.setTemperatureFromReactantsTable();
            setPressure(obj.app.mixture, pressure);

            if ~isempty(obj.app.mixture.stoichiometricMoles)
                obj.app.mixture.ratioOxidizer = obj.app.mixture.ratioOxidizer / ...
                    obj.app.mixture.stoichiometricMoles * obj.app.mixture.equivalenceRatio;
            end

            obj.updateReactantsTable();
            obj.updatePhiControls();
            obj.updateFrozenProducts();
        end

        function onEquivalenceRatioChanged(obj, event) %#ok<INUSD>
            % Rebuild reactant composition for a new equivalence ratio
            if strcmp(obj.componentValue('edit_phi', '-'), '-')
                return
            end

            if isempty(obj.componentData('UITable_R', {}))
                obj.warnInitialMixtureRequired();
                return
            end

            pressure = obj.parseComponentValue('PR2', 'first');
            equivalenceRatio = obj.parseComponentValue('edit_phi', 'first');
            previousMixture = obj.app.mixture.copy();
            obj.rebuildChemicalSystem(obj.equivalenceProductSpecies());
            obj.app.mixture = combustiontoolbox.core.Mixture(obj.app.chemicalSystem);
            obj.restoreMixtureComposition(previousMixture);
            obj.setTemperatureFromReactantsTable();
            obj.app.mixture = setProperties(obj.app.mixture, 'pressure', pressure, ...
                'equivalenceRatio', equivalenceRatio);
            obj.updateReactantsTable();
            obj.updatePhiControlsFromInput();
            obj.updateCompleteReactionProducts(equivalenceRatio);
        end

        function onReactantsTemperatureChanged(obj)
            % Apply PR1 temperature to all reactants table entries
            data = obj.componentData('UITable_R', {});

            if isempty(data)
                return
            end

            listSpecies = data(:, 1);
            temperature = {obj.parseComponentValue('PR1', 'first')};
            [temperature, FLAG_FIXED] = obj.fixSpeciesTemperatures(listSpecies, temperature, numel(listSpecies));

            if FLAG_FIXED
                data(:, 5) = temperature;
                obj.setData('UITable_R', data);
                obj.warnFixedSpeciesTemperature();
                return
            end

            data(:, 5) = repmat(temperature(1), [1, numel(listSpecies)]);
            obj.setData('UITable_R', data);
        end

        function onFrozenChemistryChanged(obj)
            % Update frozen chemistry model and product list
            if obj.hasEquilibriumSolver()
                if obj.componentValue('FrozenchemistryCheckBox', false)
                    obj.session.equilibriumSolver.caloricGasModel = ...
                        obj.session.equilibriumSolver.caloricGasModel.setThermallyPerfect();
                else
                    obj.session.equilibriumSolver.caloricGasModel = ...
                        obj.session.equilibriumSolver.caloricGasModel.setImperfect();
                end
            end

            obj.updateFrozenProducts();
        end

        function onIonizedSpeciesChanged(obj)
            % Update product list after ionized species flag changes
            obj.setValue('Products', []);
            obj.setProductItems(obj.findProductsFromReactants());
            obj.updateProductDisplayLists();
        end

        function onProductSpeciesAdded(obj)
            % Add selected database species to the products list
            items = obj.addToList(obj.componentValue('listbox_LS_DB', {}), ...
                obj.componentItems('listbox_LS'));
            obj.setItems('listbox_LS', items);
            obj.setProductItems(items);
            obj.onProductsChanged();
        end

        function onProductSpeciesRemoved(obj)
            % Remove selected species from the products list
            items = obj.removeFromList(obj.componentValue('listbox_LS', {}), ...
                obj.componentItems('listbox_LS'));
            obj.setItems('listbox_LS', items);
            obj.setProductItems(items);
            obj.onProductsChanged();
        end

        function onDatabaseSpeciesSearchChanging(obj, event)
            % Filter database species while the search field changes
            searchText = obj.eventValue(event, '');

            if isempty(searchText)
                obj.setItems('listbox_LS_DB', obj.app.database.listSpecies);
                return
            end

            obj.setItems('listbox_LS_DB', obj.filterByPrefix(obj.app.database.listSpecies, searchText));
        end

        function onDisplaySpeciesSearchChanging(obj, event)
            % Filter available and selected display species
            searchText = obj.eventValue(event, '');

            if isempty(searchText)
                obj.setItems('listbox_LS_2', obj.componentItems('listbox_LS'));
                obj.setItems('listbox_LS_display', obj.app.displaySpecies);
                return
            end

            obj.setItems('listbox_LS_2', obj.filterByPrefix(obj.componentItems('listbox_LS'), searchText));
            obj.setItems('listbox_LS_display', obj.filterByPrefix(obj.app.displaySpecies, searchText));
        end

        function onDisplaySpeciesAdded(obj)
            % Add selected species to the display list
            displayItems = obj.addToList(obj.componentValue('listbox_LS_2', {}), ...
                obj.componentItems('listbox_LS_display'));
            obj.setItems('listbox_LS_display', displayItems);
            obj.app.displaySpecies = obj.addToList(displayItems, obj.app.displaySpecies);
            obj.updateDisplaySpeciesText();
        end

        function onDisplaySpeciesRemoved(obj)
            % Remove selected species from the display list
            selectedItems = obj.componentValue('listbox_LS_display', {});
            obj.app.displaySpecies = obj.removeFromList(selectedItems, obj.app.displaySpecies);
            obj.setItems('listbox_LS_display', obj.removeFromList(selectedItems, ...
                obj.componentItems('listbox_LS_display')));
            obj.updateDisplaySpeciesText();
        end
    end

    methods (Access = private)
        function applyReactantsSelection(obj, value, event)
            preset = obj.reactantsPreset(value, event);

            if isempty(preset)
                return
            end

            obj.setValue('edit_phi', preset.phi);

            if ~isempty(preset.fuelSpecies)
                set(obj.app.mixture, preset.fuelSpecies, 'fuel', preset.fuelMoles);
            end

            if ~isempty(preset.oxidizerSpecies)
                set(obj.app.mixture, preset.oxidizerSpecies, 'oxidizer', preset.oxidizerMoles);
                obj.app.mixture.ratioOxidizer = preset.oxidizerMoles;
            end

            if ~isempty(preset.inertSpecies)
                set(obj.app.mixture, preset.inertSpecies, 'inert', preset.inertMoles);
            end
        end

        function preset = reactantsPreset(obj, value, event)
            preset = obj.emptyPreset();

            if ischar(value) || isstring(value)
                species = obj.matchSpecies(char(value));

                if isempty(species)
                    species = obj.matchSpecies(obj.eventValue(event, ''));
                end

                if isempty(species)
                    obj.showWarning('Species not found.');
                    return
                end

                obj.addCustomReactant(species, event);
                preset = [];
                return
            end

            switch value
                case 2
                    [preset.oxidizerSpecies, preset.oxidizerMoles] = obj.air();
                    preset.phi = '-';
                case 3
                    preset = obj.fuelAirPreset({'CH4'}, 1);
                case 4
                    preset = obj.fuelAirPreset({'C2H6'}, 1);
                case 5
                    preset = obj.fuelAirPreset({'C3H8'}, 1);
                case 6
                    preset = obj.fuelAirPreset({'C2H2_acetylene'}, 1);
                case 7
                    preset = obj.fuelAirPreset({'C2H4'}, 1);
                case 8
                    preset = obj.fuelAirPreset({'C6H6'}, 1);
                case 9
                    preset = obj.fuelAirPreset({'C8H18_isooctane'}, 1);
                case 10
                    preset = obj.fuelAirPreset({'H2'}, 1);
                case 11
                    preset = obj.fuelAirPreset({'CO'}, 1);
                case 12
                    preset = obj.fuelAirPreset({'CH3OH'}, 1);
                case 13
                    preset = obj.fuelAirPreset({'C2H5OH'}, 1);
                case 14
                    preset = obj.fuelAirPreset({'CH4', 'C2H6', 'C3H8'}, [0.85, 0.1, 0.05]);
                case 15
                    preset = obj.fuelAirPreset({'CO', 'H2'}, [0.5, 0.5]);
                case 16
                    preset.fuelSpecies = {'H2bLb'};
                    preset.fuelMoles = 1;
                    preset.oxidizerSpecies = {'O2bLb'};
                    preset.oxidizerMoles = 1;
                    preset.phi = '1';
                case 17
                    preset.fuelSpecies = {'RP_1'};
                    preset.fuelMoles = 1;
                    preset.oxidizerSpecies = {'O2bLb'};
                    preset.oxidizerMoles = 1;
                    preset.phi = '1';
                otherwise
                    species = obj.matchSpeciesFromEvent(event);
                    obj.addCustomReactant(species, event);
                    preset = [];
            end
        end

        function preset = emptyPreset(obj) %#ok<MANU>
            preset = struct( ...
                'fuelSpecies', {{}}, ...
                'fuelMoles', [], ...
                'oxidizerSpecies', {{}}, ...
                'oxidizerMoles', [], ...
                'inertSpecies', {{}}, ...
                'inertMoles', [], ...
                'phi', '-');
        end

        function preset = fuelAirPreset(obj, fuelSpecies, fuelMoles)
            preset = obj.emptyPreset();
            [preset.oxidizerSpecies, preset.oxidizerMoles] = obj.air();
            preset.fuelSpecies = fuelSpecies;
            preset.fuelMoles = fuelMoles;
            preset.phi = '1';
        end

        function [species, moles] = air(obj)
            [species, moles] = combustiontoolbox.utils.getAir(obj.componentValue('IdealAirCheckBox', false));
        end

        function addCustomReactant(obj, species, event) %#ok<INUSD>
            if isempty(species)
                return
            end

            if ~isempty(obj.componentData('UITable_R', {}))
                obj.setMixtureFromReactantsTable();
            end

            if any(combustiontoolbox.utils.findIndex(obj.app.mixture.listSpecies, species))
                return
            end

            set(obj.app.mixture, {species}, 'inert', 1);
        end

        function species = matchSpeciesFromEvent(obj, event)
            species = '';

            if nargin < 2 || isempty(event)
                obj.showWarning('Species not found.');
                return
            end

            species = obj.matchSpecies(obj.eventValue(event, ''));

            if isempty(species)
                obj.showWarning('Species not found.');
            end
        end

        function species = matchSpecies(obj, value)
            species = '';

            if iscell(value)
                if isempty(value)
                    return
                end

                value = value{1};
            end

            if ~(ischar(value) || isstring(value))
                return
            end

            listSpecies = obj.app.database.listSpecies;
            index = find(strcmp(listSpecies, value), 1);

            if ~isempty(index)
                species = listSpecies{index};
            end
        end

        function applyMixturePropertiesFromInputs(obj)
            temperature = obj.parseComponentValue('PR1', 'first');
            pressure = obj.parseComponentValue('PR2', 'first');

            if strcmp(obj.componentValue('edit_phi', '-'), '-')
                obj.app.mixture = setProperties(obj.app.mixture, 'temperature', temperature, 'pressure', pressure);
                return
            end

            equivalenceRatio = obj.parseComponentValue('edit_phi', 'first');
            obj.app.mixture = setProperties(obj.app.mixture, 'temperature', temperature, ...
                'pressure', pressure, 'equivalenceRatio', equivalenceRatio);
        end

        function setMixtureFromReactantsTable(obj)
            data = obj.componentData('UITable_R', {});

            if isempty(data)
                return
            end

            listSpecies = data(:, 1)';
            moles = obj.cellColumnToRowVector(data(:, 2));
            typeSpecies = data(:, 4);
            indexFuel = contains(typeSpecies, 'Fuel');
            indexOxidizer = contains(typeSpecies, 'Oxidizer');
            indexInert = contains(typeSpecies, 'Inert');

            if any(indexFuel)
                set(obj.app.mixture, listSpecies(indexFuel), 'fuel', moles(indexFuel));
            end

            if any(indexOxidizer)
                set(obj.app.mixture, listSpecies(indexOxidizer), 'oxidizer', moles(indexOxidizer));
            end

            if any(indexInert)
                set(obj.app.mixture, listSpecies(indexInert), 'inert', moles(indexInert));
            end
        end

        function restoreMixtureComposition(obj, previousMixture)
            if ~isempty(previousMixture.listSpeciesFuel)
                set(obj.app.mixture, previousMixture.listSpeciesFuel, 'fuel', previousMixture.molesFuel);
            end

            if ~isempty(previousMixture.listSpeciesOxidizer)
                set(obj.app.mixture, previousMixture.listSpeciesOxidizer, 'oxidizer', previousMixture.ratioOxidizer);
            end

            if ~isempty(previousMixture.listSpeciesInert)
                set(obj.app.mixture, previousMixture.listSpeciesInert, 'inert', previousMixture.molesInert);
            end
        end

        function setTemperatureFromReactantsTable(obj)
            data = obj.componentData('UITable_R', {});

            if isempty(data)
                return
            end

            obj.app.mixture.problemType = obj.componentValue('ProblemType', 'TP');
            speciesTemperatures = obj.cellColumnToRowVector(data(:, 5));
            [~, index] = ismember(obj.app.mixture.listSpecies, data(:, 1));
            speciesTemperatures = speciesTemperatures(index);
            setTemperatureSpecies(obj.app.mixture, speciesTemperatures);
            obj.setValue('PR1', sprintf('%.4g', obj.app.mixture.T));

            if checkTemperatureSpecies(obj.app.mixture)
                obj.warnFixedSpeciesTemperature();
            end
        end

        function updateReactantsTable(obj)
            listSpecies = [obj.app.mixture.listSpeciesInert, ...
                obj.app.mixture.listSpeciesOxidizer, ...
                obj.app.mixture.listSpeciesFuel];
            numSpecies = numel(listSpecies);

            if numSpecies == 0
                obj.clearReactantsTables();
                return
            end

            moles = [obj.app.mixture.molesInert, ...
                obj.app.mixture.molesOxidizer, ...
                obj.app.mixture.molesFuel];
            molarFraction = moles / sum(moles);
            typeSpecies = getTypeSpecies(obj.app.mixture);
            data = obj.componentData('UITable_R', {});

            if ~isempty(data) && numSpecies == size(data, 1)
                temperatures = data(:, 5)';
            else
                temperatures = repmat({obj.app.mixture.T}, [1, numSpecies]);
            end

            [temperatures, ~] = obj.fixSpeciesTemperatures(listSpecies, temperatures, numSpecies);
            tableData = [listSpecies; obj.vectorToCell(moles); ...
                obj.vectorToCell(molarFraction); typeSpecies; temperatures]';
            obj.setData('UITable_R', tableData);
            obj.setData('UITable_R2', tableData(:, 1:3));
        end

        function updatePhiControls(obj)
            if isempty(obj.app.mixture.equivalenceRatio)
                obj.setValue('edit_phi', '-');
                obj.setValue('edit_phi2', '-');
                obj.setValue('edit_phi3', '-');
                obj.setValue('edit_OF', 0);
                obj.setValue('edit_F', 0);
                return
            end

            value = sprintf('%.5g', round(obj.app.mixture.equivalenceRatio, 5));
            obj.setValue('edit_phi', value);
            obj.setValue('edit_phi2', value);
            obj.setValue('edit_phi3', value);
            obj.setValue('edit_OF', obj.app.mixture.oxidizerFuelMassRatio);
            obj.setValue('edit_F', obj.app.mixture.percentageFuel);
        end

        function updatePhiControlsFromInput(obj)
            value = obj.componentValue('edit_phi', '-');
            obj.setValue('edit_phi2', value);
            obj.setValue('edit_phi3', value);
            obj.setValue('edit_OF', obj.app.mixture.oxidizerFuelMassRatio);
            obj.setValue('edit_F', obj.app.mixture.percentageFuel);
        end

        function updateFrozenProducts(obj)
            if ~obj.componentValue('FrozenchemistryCheckBox', false)
                obj.updateProductSpecies();
                obj.updateProductDisplayLists();
                return
            end

            if contains(obj.componentValue('ProblemType', ''), 'ROCKET', 'IgnoreCase', true)
                return
            end

            data = obj.componentData('UITable_R', {});

            if isempty(data)
                obj.setProductItems({});
            else
                obj.setProductItems(data(:, 1));
            end

            obj.updateProductDisplayLists();
        end

        function updateProductSpecies(obj)
            productValue = obj.componentValue('Products', []);

            if isempty(productValue)
                obj.setProductItems(obj.findProductsFromReactants());
                return
            end

            if strcmpi(productValue, 'Complete Reaction')
                obj.onEquivalenceRatioChanged([]);
                return
            end

            FLAG_ADD = ~isempty(combustiontoolbox.utils.findIndex(obj.app.database.listSpecies, productValue));
            obj.app.chemicalSystem.setListSpecies(productValue);

            if FLAG_ADD
                obj.setProductItems(unique([obj.productSpecies(), obj.app.chemicalSystem.listSpecies], 'stable'));
            else
                obj.setProductItems(obj.app.chemicalSystem.listSpecies);
            end
        end

        function updateCompleteReactionProducts(obj, equivalenceRatio)
            if ~strcmpi(obj.componentValue('Products', []), 'Complete Reaction')
                return
            end

            obj.app.chemicalSystem.setListSpecies('complete', equivalenceRatio, ...
                obj.app.mixture.equivalenceRatioSoot);
            obj.setProductItems(obj.app.chemicalSystem.listSpecies);
            obj.updateProductDisplayLists();
        end

        function value = findProductsFromReactants(obj)
            if isempty(obj.app.mixture.listSpecies)
                value = {};
                return
            end

            listSpecies = [obj.app.mixture.listSpeciesFuel, obj.app.mixture.listSpeciesOxidizer];

            if isempty(listSpecies)
                value = obj.app.mixture.listSpecies;
                return
            end

            value = findProducts(obj.app.chemicalSystem, listSpecies, ...
                'flag_ion', obj.componentValue('IonizedspeciesCheckBox', false));
        end

        function value = productSpecies(obj)
            value = obj.componentItems('listbox_Products');

            if strcmpi(obj.componentValue('Products', ''), 'complete reaction')
                value = 'complete';
            end
        end

        function value = equivalenceProductSpecies(obj)
            if strcmpi(obj.componentValue('Products', ''), 'complete reaction')
                value = 'complete';
            else
                data = obj.componentData('UITable_R', {});
                value = data(:, 1)';
            end
        end

        function rebuildChemicalSystem(obj, listSpecies)
            if isempty(listSpecies)
                obj.app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(obj.app.database);
            else
                obj.app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(obj.app.database, listSpecies);
            end
        end

        function clearReactantsTables(obj)
            obj.setData('UITable_R', {});
            obj.setData('UITable_R2', {});
            obj.setData('UITable_P', {});
            obj.setValue('edit_phi', '-');
            obj.setValue('edit_phi2', '-');
            obj.setValue('edit_phi3', '-');
            obj.setValue('edit_OF', 0);
            obj.setValue('edit_F', 0);
        end

        function [temperature, FLAG_FIXED] = fixSpeciesTemperatures(obj, listSpecies, temperature, numSpecies)
            FLAG_FIXED = false;

            for i = 1:numSpecies
                species = obj.app.database.species.(listSpecies{i});

                if isscalar(species.T)
                    temperature{i} = species.T;
                    FLAG_FIXED = true;
                end
            end

            if FLAG_FIXED
                obj.setValue('PR1', sprintf('%.4g', max(cell2mat(temperature))));
            end
        end

        function updateProductDisplayLists(obj)
            items = obj.componentItems('listbox_Products');
            obj.setItems('listbox_LS', items);
            obj.setTitle('ListofSpeciesPanel', sprintf('List of Species - %d', obj.numProductSpecies()));
            obj.setItems('listbox_LS_2', items);
            obj.setItems('listbox_LS_display', items);
            obj.app.displaySpecies = items;
            obj.setText('text_LS', sprintf('List of Species - %d', obj.numProductSpecies()));
            obj.setText('text_LS_2', sprintf('List of Species - %d', obj.numProductSpecies()));
            obj.setText('text_LS_display', sprintf('Display Species - %d', obj.numDisplaySpecies()));
        end

        function value = numProductSpecies(obj)
            value = numel(obj.componentItems('listbox_Products'));
        end

        function value = numDisplaySpecies(obj)
            value = numel(obj.componentItems('listbox_LS_display'));
        end

        function updateDisplaySpeciesText(obj)
            obj.setText('text_LS_display', sprintf('Display Species - %d', obj.numDisplaySpecies()));
        end

        function values = filterByPrefix(obj, items, searchText) %#ok<INUSD>
            values = {};

            if isempty(items)
                return
            end

            searchText = char(searchText);
            index = false(1, numel(items));

            for i = 1:numel(items)
                index(i) = startsWith(items{i}, searchText, 'IgnoreCase', false);
            end

            values = items(index);
        end

        function items = addToList(obj, values, items)
            values = obj.asCell(values);
            items = obj.asCell(items);

            if isempty(values)
                return
            end

            items = unique([items, values], 'stable');
        end

        function items = removeFromList(obj, values, items)
            values = obj.asCell(values);
            items = obj.asCell(items);

            if isempty(values)
                return
            end

            items = items(~ismember(items, values));
        end

        function values = asCell(obj, values) %#ok<INUSD>
            if isempty(values)
                values = {};
            elseif ischar(values)
                values = {values};
            elseif isstring(values)
                values = cellstr(values(:))';
            elseif ~iscell(values)
                values = {values};
            else
                values = values(:)';
            end
        end

        function value = parseComponentValue(obj, componentName, varargin)
            value = combustiontoolbox.gui.AppInput.parseValue(obj.componentValue(componentName, []), varargin{:});
        end

        function value = componentValue(obj, componentName, defaultValue)
            value = defaultValue;

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Value')
                value = obj.app.(componentName).Value;
            end
        end

        function value = eventValue(obj, event, defaultValue) %#ok<INUSD>
            value = defaultValue;

            if isempty(event)
                return
            end

            if isstruct(event) && isfield(event, 'Value')
                value = event.Value;
                return
            end

            if isobject(event) && isprop(event, 'Value')
                value = event.Value;
            end
        end

        function value = componentItems(obj, componentName)
            value = {};

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Items')
                value = obj.app.(componentName).Items;
            end
        end

        function data = componentData(obj, componentName, defaultValue)
            data = defaultValue;

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Data')
                data = obj.app.(componentName).Data;
            end
        end

        function setValue(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Value')
                obj.app.(componentName).Value = value;
            end
        end

        function setItems(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Items')
                value = obj.asCell(value);
                obj.app.(componentName).Items = value;
            end
        end

        function setProductItems(obj, value)
            obj.setItems('listbox_Products', value);
        end

        function setData(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Data')
                obj.app.(componentName).Data = value;
            end
        end

        function setText(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Text')
                obj.app.(componentName).Text = value;
            end
        end

        function setTitle(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Title')
                obj.app.(componentName).Title = value;
            end
        end

        function value = hasComponent(obj, componentName)
            value = isobject(obj.app) && isprop(obj.app, componentName) && ~isempty(obj.app.(componentName));
        end

        function value = hasEquilibriumSolver(obj)
            value = ~isempty(obj.session) && isprop(obj.session, 'equilibriumSolver') ...
                && ~isempty(obj.session.equilibriumSolver);
        end

        function value = cellColumnToRowVector(obj, values) %#ok<INUSD>
            if iscell(values)
                value = cell2mat(values)';
            else
                value = values';
            end
        end

        function value = vectorToCell(obj, values) %#ok<INUSD>
            value = num2cell(values);
        end

        function warnInitialMixtureRequired(obj)
            obj.setWarningStatus();
            obj.showWarning('First, define initial mixture');
        end

        function warnFixedSpeciesTemperature(obj)
            obj.showWarning(['The species selected as reactants can only be evaluated at its defined temperature.' newline ...
                'The temperature shown as the temperature of the reactant is a ficticious value! The species are unmixed.']);
        end

        function showWarning(obj, message)
            if obj.shouldShowAlert()
                uialert(obj.app.UIFigure, {message}, 'Warning', 'Icon', 'warning');
            elseif obj.hasComponent('consolePanel')
                obj.app.consolePanel.setOutput(message);
            elseif obj.hasComponent('Console_text') && isprop(obj.app.Console_text, 'Value')
                obj.app.Console_text.Value = message;
            else
                warning('AppSpeciesPanel:Warning', '%s', message);
            end
        end

        function value = shouldShowAlert(obj)
            value = obj.hasComponent('UIFigure') ...
                && isvalid(obj.app.UIFigure) ...
                && strcmpi(char(obj.app.UIFigure.Visible), 'on');
        end

        function setWarningStatus(obj)
            if obj.hasComponent('statusPanel')
                obj.app.statusPanel.setWarning();
            end
        end

    end
end
