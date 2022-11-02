function index = find_ind(LS, species)
    % Find the index of the species based on the given list (LS)
    %
    % Args:
    %     LS (cell):      List of species
    %     species (cell): Species to find index values
    %
    % Returns:
    %     index (float):  List with the index of the species based on the given list (LS)

    try
        % Initialization
        index = [];
        % Check inputs
        [LS, FLAG_1] = check_iscell(LS);
        [species, FLAG_2] = check_iscell(species);
        % Find index
        if FLAG_1 && FLAG_2
            [bool_index, loc] = ismember(LS, species);
            [~, index_sort] = sort(loc(bool_index));
            index = find(bool_index);
            index = index(index_sort);
        end

    catch ME
        print_error(ME, 'type', 'Warning', 'Solution', 'Returning an empty index value.')
    end

end

% SUB-PASS FUNCTIONS
function [variable, FLAG] = check_iscell(variable)
    % Convert to cell if the given variable is not a cell
    if ~isempty(variable)
        FLAG = true;
    else
        FLAG = false;
        return
    end
    
    if ~iscell(variable)
        variable = {variable};
    end

end
