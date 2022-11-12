function gui_clear_results(app)
    % Function that clears the result tab panel, setting them to 0
    
    % Clear Tab results: parameters
    update_properties_common(app, 0, '_1');
    update_properties_common(app, 0, '_2');
    update_properties_common(app, 0, '_3');
    update_properties_common(app, 0, '_4');
    update_properties_common(app, 0, '_5');
    update_properties_rocket(app, 0, '_2');
    update_properties_rocket(app, 0, '_3');
    update_properties_rocket(app, 0, '_4');
    update_properties_rocket(app, 0, '_5');
    % Clear error
    app.text_error_problem.Value = 0;
end

function update_properties_common(app, value, suffix)
    % Update properties to the given value
    app.(['text_T', suffix]).Value = value;
    app.(['text_p', suffix]).Value = value;
    app.(['text_r', suffix]).Value = value;
    app.(['text_h', suffix]).Value = value;
    app.(['text_e', suffix]).Value = value;
    app.(['text_cp', suffix]).Value = value;
    app.(['text_s', suffix]).Value = value;
    app.(['text_gamma', suffix]).Value = value;
    app.(['text_W', suffix]).Value = value;
    app.(['text_sound', suffix]).Value = value;
    app.(['text_u', suffix]).Value = value;
    app.(['text_M', suffix]).Value = value;
end

function update_properties_rocket(app, value, suffix)
    % Update rocket propellant performance parameters
    app.(['text_Aratio', suffix]).Value = value;
    app.(['text_Cstar', suffix]).Value = value;
    app.(['text_Ivac', suffix]).Value = value;
    app.(['text_Isp', suffix]).Value = value;
end