function LS = gui_value2list(app, value, LS, action)
    % Add/remove (action) selected value to/from the list of species (LS)
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     value (cell): Cell array of char containing the selected species
    %     LS (cell): Cell array of char containing the list of species
    %     action (char): 'add' or 'remove'
    %
    % Returns:
    %     LS (cell): Cell array of char containing the list of species after the action
    %
    % Note:
    %     If value is empty, the function returns LS without any change

    if isempty(value)
        return
    end

    switch lower(action)
        case 'add'
            LS = add_value(value, LS);
        case 'remove'
            LS = remove_value(value, LS);
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