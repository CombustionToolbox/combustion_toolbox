function LS_output = gui_add_value2list(obj, value, LS_output)
    % Add selected value to LS_output

    % Check if the value is already included in the output list
    n_pass = [];
    for n = length(value):-1:1
        if ~any(strcmpi(LS_output, value{n}))
            n_pass = [n_pass, n];
        end
    end
    LS_output = [LS_output, value(n_pass)];
end