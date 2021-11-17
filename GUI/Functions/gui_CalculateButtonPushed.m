function obj = gui_CalculateButtonPushed(obj, event)
    % Solve selected problem, update GUI with the results and generate
    % predefined plots.
    try
        % Set lamp to Working color
        obj.Lamp.Color = obj.color_lamp_working;
        % Get List of Species considered (reactants + products)
        obj = get_listSpecies_gui(obj);
        % Initialize variable app
        app = App('fast', obj.DB_master, obj.DB, obj.LS);
        % Set FLAG GUI
        app.Misc.FLAG_GUI = true;
        % Get initial conditions
        app = get_input_constrains(obj, app);
        app = gui_get_reactants(obj, event, app);
        % Solve selected problem
        app = SolveProblem(app, app.PD.ProblemType);
        % Save results
        results = save_results(obj, app);
        % Update GUI with the last results of the set
        obj = update_results_gui(obj, results);
        % Display results (plots)
        postResults(app);
        % Set lamp to Done color
        obj.Lamp.Color = obj.color_lamp_done;
    catch ME
        % Set lamp to Error color
        obj.Lamp.Color = obj.color_lamp_error;
        % Print error
        errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
        fprintf('%s\n', errorMessage);
        uiwait(warndlg(errorMessage));
    end
end

% SUB-PASS FUNCTIONS
function obj = update_results_gui(obj, app)
    % Update GUI with the results computed
end

function obj = get_listSpecies_gui(obj)
    obj.LS = obj.listbox_LS.Value;
    if isempty(obj.LS)
        % Get default value
        obj.LS = 'Soot Formation';
    end
end
function results = save_results(obj, app)
    % Save results in the UITree
    % 1. Save data 
    results.mix1 = app.PS.strR;
    results.mix2 = app.PS.strP;
    results.ProblemType = obj.ProblemType.Value;
    results.reaction = obj.Reaction.Value;
    results.LS = obj.LS;
    results.LS_products = obj.LS_products;
    % 2. Save results in the UITree 
end

function app = get_input_constrains(obj, app)
    [app.PD.TR.value, app.FLAG_PR1] = gui_get_prop(obj, 'TR', obj.PR1.Value);
    [app.PD.pR.value, app.FLAG_PR2] = gui_get_prop(obj, 'pR', obj.PR2.Value);
    [app.PD.phi.value, app.FLAG_phi] = gui_get_prop(obj, 'phi', obj.edit_phi.Value);
    app.PD.ProblemType = obj.ProblemType.Value;
    switch app.PD.ProblemType
        case 'TP' % * TP: Equilibrium composition at defined T and p
            [app.PD.TP.value, app.FLAG_PP1] = gui_get_prop(obj, 'TP', obj.PP1.Value);
            [app.PD.pP.value, app.FLAG_PP2] = gui_get_prop(obj, 'pP', obj.PP2.Value);
        case 'HP' % * HP: Adiabatic T and composition at constant p
            [app.PD.pP.value, app.FLAG_PP2] = gui_get_prop(obj, 'pP', obj.PP2.Value);
        case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            [app.PD.pP.value, app.FLAG_PP2] = gui_get_prop(obj, 'pP', obj.PP2.Value);
            app.PD.phi.value = 1*ones(1, length(app.PD.pP.value));
        case 'TV' % * TV: Equilibrium composition at defined T and constant v
            [app.PD.TP.value, app.FLAG_PP1] = gui_get_prop(obj, 'TP', obj.PP1.Value);
            [app.PD.pP.value, app.FLAG_PP2] = gui_get_prop(obj, 'pP', obj.PP2.Value);
        case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
            [app.PD.pP.value, app.FLAG_PP2] = gui_get_prop(obj, 'pP', obj.PP2.Value);
        case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            app.PD.vP_vR.value = gui_get_prop(obj, 'vP_vR', obj.PP2.Value);
            app.PD.phi.value = 1*ones(1, length(app.PD.vP_vR.value));
        case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE

        case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE

        case 'DET' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            % No additional constrains
        case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION

    end
end