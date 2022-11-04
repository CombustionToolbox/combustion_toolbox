function app = gui_get_parameters(app)
    % GET EQUIVALENCE RATIO
    [app.PD.phi.value, app.flag_phi] = gui_get_prop(app, app.edit_phi.Value, 'phi');
    if isempty(app.PD.phi.value)
        app.PD.phi.value = 1;
    end
    % GET CONDITIONS
    [app.PR1_vector, app.flag_PR1] = gui_get_prop(app, app.PR1.Value, app.PR1_var_name);
    [app.PR2_vector, app.flag_PR2] = gui_get_prop(app, app.PR2.Value, app.PR2_var_name);
    [app.PP1_vector, app.flag_PP1] = gui_get_prop(app, app.PP1.Value, app.PP1_var_name);
    [app.PP2_vector, app.flag_PP2] = gui_get_prop(app, app.PP2.Value, app.PP2_var_name);
    if strcmpi(app.PD.ProblemType, 'DET_OVERDRIVEN')
        [app.PR3_vector, app.flag_PR3] = gui_get_prop(app, app.PR3.Value, app.PR3_var_name);
    end
end