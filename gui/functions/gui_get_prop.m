function [value, flag_prop] = gui_get_prop(app, name, value, varargin)
    if value(1) == '['
        flag_prop = true;
        value = sscanf(value, '[%f:%f:%f]');
        value = value(1):value(2):value(3);
    else
        flag_prop = false;
        value = sscanf(value,'%f');
    end

    extra_parameters = nargin - 3;
    if extra_parameters
        n = 1;
        direction = 'first';
        for i = 1:2:extra_parameters
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
end