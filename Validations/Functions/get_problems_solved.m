function value = get_problems_solved(varargin)
    % Get problems solved based on the length of the given variable
    mix = varargin{1};
    variable = varargin{2};
    if nargin == 3
        subvariable = varargin{3};
        for i=length(mix):-1:1
            value(i) = length(mix{i}.(variable).(subvariable));
        end
    else
        for i=length(mix):-1:1
            value(i) = length(mix{i}.(variable));
        end
    end
    value = sum(value);
end