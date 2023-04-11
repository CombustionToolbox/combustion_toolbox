function gui_compute_mach_or_velocity(app, inputname)
    % Function that computes the pre-shock Mach number of the mixture from 
    % a given pre-shock velocity or viceversa.

    % If the problem type is not a shock problem, return
    if ~any(contains(app.ProblemType.Value, 'SHOCK', 'IgnoreCase', true))
        return
    end

    % If the table is empty, return a dash
    if isempty(app.UITable_R.Data)
        app.PR3.Value = '-';
        return
    end

    % Compute the pre-shock Mach number/velocity
    temp_app = gui_create_temp_app(app, [], false);
    mix1 = temp_app.PS.strR{1};

    if strcmpi(inputname, 'Mach')
        [M1, FLAG] = gui_get_prop(app, 'M1', app.PR4.Value);
        u1 = M1 * soundspeed(mix1);
        app.PR3.Value = compute_vector_or_scalar(u1, FLAG);
    else
        [u1, FLAG] = gui_get_prop(app, 'u1', app.PR3.Value);
        M1 = u1 / soundspeed(mix1);
        app.PR4.Value = compute_vector_or_scalar(M1, FLAG);
    end

end

% SUB-PASS FUNCTIONS
function value_char = compute_vector_or_scalar(value, FLAG)
    if FLAG
        value_char = sprintf('%.2f:%.2f:%.2f', value(1), value(2) - value(1), value(end));
    else
        value_char = sprintf('%.2f', value);
    end

end