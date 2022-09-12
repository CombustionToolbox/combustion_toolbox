function index = find_ind(LS, species)
    % Find the index of the species based on the given list (LS)
    %
    % Args:
    %     LS (cell):      List of species
    %     species (cell): Species to find index values
    %
    % Returns:
    %     index (float):  List with the index of the species based on the given list (LS)
    
    % Check inputs
    check_iscell(LS);
    check_iscell(species);
    % Find index
    [bool_index, loc] = ismember(LS, species);
    [~, index_sort] = sort(loc(bool_index));
    index = find(bool_index);
    index = index(index_sort);
end

% SUB-PASS FUNCTIONS
function variable = check_iscell(variable)
    % Convert to cell if the given variable is not a cell
    if ~iscell(variable)
        variable = {variable};
    end
end