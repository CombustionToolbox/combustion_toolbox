function listSpecies = findSpecies(listSpecies, cond_with, type_with, cond_without, type_without)
    % Find species in the given list that contain all/any
    % elements of cond_with and that not include all/any elements of
    % cond_without
    %
    % Args:
    %     listSpecies (cell): List of species
    %     cond_with (cell): List of elements to include
    %     type_with (char): Satisfy all or any of the elements in cond_with
    %     cond_without (cell): List of elements to avoid
    %     type_without (char): Satisfy all or any of the elements in cond_without
    %
    % Returns:
    %     listSpecies (cell): List of species
    %
    % Examples:
    %     * listSpecies = findSpecies(listSpecies, {'C','N','O','minus','plus','Ar'}, 'any',...
    %                                {'I', 'S', 'L', 'T', 'P', 'F', 'ab', 'W',...
    %                                 'Z','X','R','Os','Cr','H','Br','G','K',...
    %                                 'U','Co','Cu','B','V','Ni','Na','Mg',...
    %                                 'Mo','Ag','Nb','Cb','Cl','D','T',...
    %                                 'Ca','Cs','Ne','Cd','Mn'}, 'all')
    %     * listSpecies = findSpecies(DB.listSpecies, {}, 'any', {'_M'}, 'all')

    % Initialization
    FLAG_WITH = [];
    FLAG_WITHOUT = [];
    % Get all flags
    for i = length(cond_with):-1:1
        FLAG_WITH(:, i) = contains(listSpecies, cond_with{i});
    end

    for i = length(cond_without):-1:1
        FLAG_WITHOUT(:, i) = ~contains(listSpecies, cond_without{i});
    end

    % Get functions all or any for each list of elements
    f_with = get_fun_type(type_with);
    f_without = get_fun_type(type_without);
    % Obtain flags that contain (not contain) all/any of the given elements
    index = get_FLAG(f_with, FLAG_WITH) & get_FLAG(f_without, FLAG_WITHOUT);
    % Get list of species
    listSpecies = listSpecies(index);
end

% SUB-PASS FUNCTION
function f = get_fun_type(type)
    % Get function all or any
    if isempty(type)
        f = [];
        return
    end

    switch lower(type)
        case 'all'
            f = @all;
        case 'any'
            f = @any;
        otherwise
            f = [];
    end

end

function FLAG = get_FLAG(f, FLAG)
    % Get flag that satisfy all/any of the given elements
    if ~isempty(FLAG)
        FLAG = f(FLAG, 2);
    else
        FLAG = true;
    end

end
