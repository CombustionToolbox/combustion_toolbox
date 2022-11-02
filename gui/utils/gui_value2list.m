function LS = gui_value2list(obj, value, LS, action)
    % Add/remove (action) selected value to LS
    if ~isempty(value)
        switch lower(action)
            case 'add'
                LS = add_value(value, LS);
            case 'remove'
                LS = remove_value(value, LS);
        end
    end
end

% SUB-PASS FUNCTIONS
function LS = add_value(value, LS)
    % Add value to LS

    % Check if the value is already included in the list
    n_pass = [];
    for n = length(value):-1:1
        if ~any(strcmp(LS, value{n}))
            n_pass = [n_pass, n];
        end
    end
    LS = [LS, value(n_pass)];
end

function LS = remove_value(value, LS)
    % Remove value from LS

    % Remove the species from the list
    N = length(LS);
    n_pass = true(1, N);
    for n = N:-1:1
        if any(strcmp(value, LS{n}))
            n_pass(n) = false;
        end
    end
    LS = LS(n_pass);
end