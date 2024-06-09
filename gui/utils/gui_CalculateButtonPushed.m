function app = gui_CalculateButtonPushed(app, event)
    % Solve selected problem, update GUI with the results, and generate
    % predefined plots
    %
    % Args:
    %   app (object): Combustion Toolbox app object
    %   event (object): Event object
    %
    % Returns:
    %   app (object): Combustion Toolbox app object
    
    % Definitions
    propertyNameR1 = 'temperature';
    propertyNameR2 = 'pressure';
    equivalenceRatio = gui_get_prop(app.edit_phi.Value);
    problemType = app.ProblemType.Value;
    listSpecies = app.listbox_Products.Items;
    FLAG_EQUIVALENCE_RATIO = ~isempty(equivalenceRatio);
    FLAG_RESULTS = app.PrintresultsCheckBox.Value;

    try
        % Set lamp to Working color
        app.Lamp.Color = app.color_lamp_working;
        
        % Get properties
        propertyR1 = gui_get_prop(app.PR1.Value);
        propertyR2 = gui_get_prop(app.PR2.Value);

        switch problemType
            case {'TP', 'TV'}
                propertyP1 = gui_get_prop(app.PP1.Value);
                propertyP2 = gui_get_prop(app.PP2.Value);
            otherwise
                propertyP1 = propertyR1;
                propertyP2 = propertyR2;
        end
        

        % Check problem type
        if strcmpi(app.ProblemType.Value(2), 'V')
            propertyNameR2 = 'volume';
        end
        
        % Define chemical system
        app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database, listSpecies);

        % Temporal mixture
        tempMixture = app.mixture.copy();
        
        % Initialize mixture
        app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);

        % Define chemical state
        if ~isempty(tempMixture.listSpeciesFuel)
            set(app.mixture, tempMixture.listSpeciesFuel, 'fuel', tempMixture.molesFuel);
        end

        if ~isempty(tempMixture.listSpeciesOxidizer)
            set(app.mixture, tempMixture.listSpeciesOxidizer, 'oxidizer', tempMixture.ratioOxidizer);
        end

        if ~isempty(tempMixture.listSpeciesInert)
            set(app.mixture, tempMixture.listSpeciesInert, 'inert', tempMixture.molesInert);
        end
        
        % Set ratioOxidizer
        app.mixture.ratioOxidizer = tempMixture.ratioOxidizer;
        
        % Define mixArray
        if FLAG_EQUIVALENCE_RATIO
            mixArray1 = setProperties(app.mixture, propertyNameR1, propertyR1, propertyNameR2, propertyR2, 'equivalenceRatio', equivalenceRatio);
            mixArray2 = setProperties(app.mixture, propertyNameR1, propertyP1, propertyNameR2, propertyP2, 'equivalenceRatio', equivalenceRatio);
        else
            mixArray1 = setProperties(app.mixture, propertyNameR1, propertyR1, propertyNameR2, propertyR2);
            mixArray2 = setProperties(app.mixture, propertyNameR1, propertyP1, propertyNameR2, propertyP2);
        end

        % Update GUI terminal
        gui_update_terminal(app, 'start');

        % Select solver and solve problem
        switch problemType
            case {'TP', 'HP', 'SP', 'TV', 'EV', 'SV'}
                % Select solver
                solver = combustiontoolbox.equilibrium.EquilibriumSolver('problemType', problemType, 'FLAG_RESULTS', FLAG_RESULTS);
                % Solve problem
                solver.solveArray(mixArray2);
            case {'SHOCK_I'}
                solver = combustiontoolbox.shockdetonation.ShockSolver('problemType', problemType, 'FLAG_RESULTS', FLAG_RESULTS);
            case {'DET'}
                solver = combustiontoolbox.shockdetonation.DetonationSolver('problemType', problemType, 'FLAG_RESULTS', FLAG_RESULTS);
            case {'ROCKET'}
                solver = combustiontoolbox.shockdetonation.RocketSolver('problemType', problemType, 'FLAG_RESULTS', FLAG_RESULTS);
            otherwise
                error('Problem type %s is not found', problemType);
        end

        % Update GUI terminal
        gui_update_terminal(app, 'finish', solver.time);

        % Save results
        [results, app.temp_results] = save_results(app, problemType, mixArray1, mixArray2);
        
        % Update GUI with the last results of the set
        gui_update_results(app, results);
        
        % Display results (plots)
        switch lower(app.Report_type.Value)
            case {'auto'}
                solver.plot(mixArray2);
        end

        % % Get List of Species considered (reactants + products)
        % app = get_listSpecies_gui(app);
        % 
        % % Initialize variable self
        % self = App('fast', app.DB_master, app.DB, app.LS);
        % 
        % % Get Miscellaneous
        % self = get_miscellaneous_values(app, self);
        % 
        % % Get tuning values and constants
        % self = get_tuning_values(app, self);
        % 
        % % Get initial conditions
        % self = get_input_constrains(app, self);
        % self = gui_get_reactants(app, event, self);
        % 
        % % Get display species
        % self = gui_get_display_species(app, self);
        % 
        % % Update GUI terminal
        % gui_update_terminal(app, self, 'start');
        % 
        % % Solve selected problem
        % self = solve_problem(self, self.PD.ProblemType);
        % 
        % % Save results
        % [results, app.temp_results] = save_results(app, self);
        % 
        % % Update GUI with the last results of the set
        % gui_update_results(app, results);
        % 
        % % Update GUI custom figures tab
        % gui_update_custom_figures(app);
        % 
        % % Update GUI terminal
        % gui_update_terminal(app, self, 'finish')
        % 
        % % Display results (plots)
        % switch lower(app.Report_type.Value)
        %     case {'auto'}
        %         post_results(self);
        % end

        % Set lamp to Done color
        app.Lamp.Color = app.color_lamp_done;

    catch ME
        % Set lamp to Error color
        app.Lamp.Color = app.color_lamp_error;
        % Print error
        message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message)};
        uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
    end

end

% SUB-PASS FUNCTIONS


function self = get_miscellaneous_values(obj, self)
    % Get miscellaneous properties from GUI and set into self
    self.Misc = obj.Misc;
end

function self = get_tuning_values(obj, self)
    % Get properties from GUI and set into self
    self.PD = obj.PD;
    self.TN = obj.TN;
    self.C = obj.C;
end

function gui_update_results(obj, results)
    % Update:
    %  1. GUI with the last results computed
    %  2. GUI-UITree with all the results
    
    % Update GUI with the last result computed
    gui_write_results(obj, results, 1);
    
    % Update UITree with all the results
    gui_add_nodes(obj.Node_Results, results)
end

function obj = get_listSpecies_gui(obj)
    obj.LS = obj.listbox_Products.Items;
    if isempty(obj.LS)
        % Get default value
        obj.LS = list_species('Soot Formation');
    end

end

function [results, temp_results] = save_results(obj, problemType, mixArray1, mixArray2, varargin)
    % Save results
    
    % Definitions
    numCases1 = length(mixArray1);
    numCases2 = length(mixArray2);
    numCases = max(numCases1, numCases2);

    % % 
    % if strcmpi(self.PD.ProblemType, 'ROCKET')
    %     label_mix1 = 'strR';
    %     if self.PD.FLAG_IAC
    %         label_mix2 = 'mix2_c';
    %         label_mix3 = 'mix3';
    %         label_mix4 = 'strP';
    %         label_mix5 = [];
    %     else
    %         label_mix2 = 'str2';
    %         label_mix3 = 'mix2_c';
    %         label_mix4 = 'mix3';
    %         label_mix5 = 'strP';
    %     end
    % 
    % else
    %     label_mix1 = 'strR';
    %     label_mix2 = 'strP';
    % end

    for i = numCases:-1:1
        
        if numCases1 > numCases2
            results(i).mix1 = mixArray1(i);
            results(i).mix2 = mixArray2;
        elseif numCases1 < numCases2
            results(i).mix1 = mixArray1;
            results(i).mix2 = mixArray2(i);
        else
            results(i).mix1 = mixArray1(i);
            results(i).mix2 = mixArray2(i);
        end

        results(i).ProblemType = problemType;

        if ischar(obj.Reactants.Value)
            results(i).Reactants = 'Custom';
        else
            results(i).Reactants = obj.Reactants.Items{obj.Reactants.Value};
        end

        results(i).Products = obj.Products.Value;

        if isempty(results(i).Products)
            results(i).Products = 'Default';   
        end

        results(i).listSpecies = mixArray2(i).chemicalSystem.listSpecies;
        results(i).listProducts = mixArray2(i).chemicalSystem.listProducts;
        results(i).UITable_R_Data = obj.UITable_R.Data;
    end

    % if strcmpi(self.PD.ProblemType, 'ROCKET')
    %     try
    %         label_mix = {label_mix3, label_mix4, label_mix5};
    %         for mix = label_mix
    %             for i = numCases:-1:1
    %                 results(i).(mix{:}) = self.PS.(mix{:}){i};
    %             end
    % 
    %         end
    % 
    %     catch
    %         % Nothing to do
    %     end
    % 
    % end

    % Save temporally this parametric study
    temp_results = results(i);
end

function self = get_input_constrains(obj, self)
    % Get input constrains
    self.PD.ProblemType = obj.ProblemType.Value;
    self = set_prop(self, 'TR', obj.PR1.Value, 'pR', obj.PR2.Value, 'phi', obj.edit_phi.Value);
    switch self.PD.ProblemType
        case 'TP' % * TP: Equilibrium composition at defined T and p
            self = set_prop(self, 'TP', obj.PP1.Value, 'pP', obj.PP2.Value);
        case 'HP' % * HP: Adiabatic T and composition at constant p
            self = set_prop(self, 'pP', obj.PP2.Value);
        case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            self = set_prop(self, 'pP', obj.PP2.Value);
        case 'TV' % * TV: Equilibrium composition at defined T and constant v
            self = set_prop(self, 'TP', obj.PP1.Value, 'pP', obj.PP2.Value);
        case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
            self = set_prop(self, 'pP', obj.PP2.Value);
        case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            self = set_prop(self, 'vP_vR', obj.PP4.Value);
        case {'SHOCK_I', 'SHOCK_R', 'SHOCK_POLAR'} % * SHOCK_I, SHOCK_R, and SHOCK_POLAR
            self = set_prop(self, 'u1', obj.PR3.Value);
        case {'SHOCK_OBLIQUE', 'SHOCK_OBLIQUE_R', 'SHOCK_POLAR_R'} % * SHOCK_OBLIQUE, SHOCK_OBLIQUE_R, and SHOCK_POLAR
            self = set_prop(self, 'u1', obj.PR3.Value);

            if ~isempty(obj.PR5.Value)
                self = set_prop(self, 'beta', obj.PR5.Value);
            else
                self = set_prop(self, 'theta', obj.PP5.Value);
            end

        case {'DET', 'DET_R'} % * DET and DET_R
            % No additional constrains
        case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R', 'DET_POLAR'} % * DET_OVERDRIVEN, DET_UNERDRIVEN, and DET_POLAR
            self = set_prop(self, 'drive_factor', obj.PR3.Value);
        case {'DET_OBLIQUE', 'DET_OBLIQUE_R'} % * DET_OBLIQUE and DET_OBLIQUE_R
            self = set_prop(self, 'drive_factor', obj.PR3.Value);
            
            if ~isempty(obj.PR4.Value)
                self = set_prop(self, 'beta', obj.PR4.Value);
            else
                self = set_prop(self, 'theta', obj.PP4.Value);
            end

        case {'ROCKET'} % * ROCKET: ROCKET PROPELLANT PERFORMANCE
            % Get model
            self.PD.FLAG_IAC = obj.FLAG_IAC.Value;
            if ~obj.FLAG_IAC.Value
                if ~isempty(obj.PP1.Value)
                    self = set_prop(self, 'Aratio_c', obj.PP1.Value);
                elseif ~isempty(obj.PP2.Value)
                    self = set_prop(self, 'mass_flux', obj.PP2.Value);
                else
                    % Set lamp to Error color
                    obj.Lamp.Color = obj.color_lamp_error;
                    % Print error
                    message = {'The FAC model needs an additional value! The contraction factor A_chamber/A_throat or the mass flux.'};
                    uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
                end

            end

            if ~isempty(obj.PR3.Value)
                % Set Aratio combustor to throat region (subsonic region)
                self = set_prop(self, 'Aratio', obj.PR3.Value);
                self.PD.FLAG_SUBSONIC = true;
            end

            if ~isempty(obj.PP3.Value)
                % Set Aratio throat to exit region (supersonic region)
                self = set_prop(self, 'Aratio', obj.PP3.Value);
                self.PD.FLAG_SUBSONIC = false;
            end

    end
    if contains(obj.Products.Value, 'complete', 'IgnoreCase', true)
        self.S.FLAG_COMPLETE = true;
    end

end

function gui_update_custom_figures(obj)
    % Update GUI custom figures tab
    try
        % Remove previous nodes from UITree
        delete(obj.Mixtures.Children);
        delete(obj.Variable_x.Children);
        delete(obj.Variable_y.Children);
        % Add new nodes
        results = obj.temp_results;
        fields = fieldnames(results);
        FLAG_fields = contains(fields, 'str') | contains(fields, 'mix');
        fields = fields(FLAG_fields);
        % Mixtures
        add_node(obj, 'Mixtures', fields);
        % Variables x
        try
            fields = sort(fieldnames(obj.temp_results.strP{1}));
        catch
            fields = sort(fieldnames(obj.temp_results.mix3{1}));
        end

        add_node(obj, 'Variable_x', fields);
        % Variables y
        add_node(obj, 'Variable_y', fields);

    catch ME
        % Print error
        fprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
    end
    
end


function add_node(obj, field_master, field_names)
    % Add node to the UItree
    for i = 1:length(field_names)
        uitreenode(obj.(field_master), 'Text', field_names{i});
    end

end

function self = gui_get_display_species(obj, self)
    % Set display species in variable self from the display species itembox
    % (GUI)
    self.Misc.display_species = obj.listbox_LS_display.Items;
end