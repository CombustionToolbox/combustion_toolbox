function obj = gui_CalculateButtonPushed(obj)
    % Solve selected problem, update GUI with the results and generate
    % predefined plots.
    try
        % Initialize
        app = App('Soot formation');
        % Get initial conditions
        app = get_input_constrains(obj, app);
        app = get_current_reactants_gui(obj, app);
        % Solve selected problem
        app = SolveProblem(app, app.PD.ProblemType);
        % Display results (plots)
        postResults(app);
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end

% SUB-PASS FUNCTIONS
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