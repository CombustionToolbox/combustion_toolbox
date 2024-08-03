function listSpecies = gui_value2list(app, value, listSpecies, action)
    % Add/remove (action) selected value to/from the list of species
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     value (cell): Cell array of char containing the selected species
    %     listSpecies (cell): Cell array of char containing the list of species
    %     action (char): 'add' or 'remove'
    %
    % Returns:
    %     listSpecies (cell): Cell array of char containing the list of species after the action
    %
    % Note:
    %     If value is empty, the function returns listSpecies without any change

    if isempty(value)
        return
    end

    switch lower(action)
        case 'add'
            listSpecies = add_value(value, listSpecies);
        case 'remove'
            listSpecies = remove_value(value, listSpecies);
    end

end

% SUB-PASS FUNCTIONS
function listSpecies = add_value(value, listSpecies)
    % Add value to LS
    
    % Definitions
    n_pass = [];

    % Check if the value is already included in the list
    for n = length(value):-1:1

        if ~any(strcmp(listSpecies, value{n}))
            n_pass = [n_pass, n];
        end

    end

    listSpecies = [listSpecies, value(n_pass)];
end

function listSpecies = remove_value(value, listSpecies)
    % Remove value from LS
    
    % Definitions
    numSpecies = length(listSpecies);
    n_pass = true(1, numSpecies);

    % Remove the species from the list
    for n = numSpecies:-1:1
        
        if any(strcmp(value, listSpecies{n}))
            n_pass(n) = false;
        end

    end
    
    listSpecies = listSpecies(n_pass);
end