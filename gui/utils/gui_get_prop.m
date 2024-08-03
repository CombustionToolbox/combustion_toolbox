function [value, FLAG_ARRAY] = gui_get_prop(value, varargin)
    % GUI routine to get a property value from a char
    %
    % Args:
    %     value (char): Property value as a char
    %
    % Optional Args:
    %     * number (float): number of values to get
    %     * direction (char): 'first' or 'last' to get the first or last value
    %
    % Returns:
    %     value (float): Property value as a float
    %     FLAG_ARRAY (bool): True if the value is an array
    
    % Initialization
    FLAG_ARRAY = false;

    % Check inputs
    if isempty(value)
        return
    end
    
    % Read value
    if value(1) == '['
        FLAG_ARRAY = true;
        value = str2num(value);
    elseif sum(value == ':') == 2
        FLAG_ARRAY = true;
        value = sscanf(value, '%f:%f:%f');
        value = value(1):value(2):value(3);
    else
        FLAG_ARRAY = false;
        value = sscanf(value, '%f');
    end

    if nargin == 1
        return
    end

    n = 1;
    direction = 'first';
    for i = 1:2:nargin-1

        switch varargin{i}
            case {'number', 'n', 'nindex'}
                n = varargin{i + 1};
            case {'direction', 'get', 'position', 'pos'}
                direction = varargin{i + 1};
        end

    end

    index = find(value, n, direction);
    value = value(index);
end