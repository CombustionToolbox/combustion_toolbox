function C = vector2cell(value)
    % Create cell array from vector
    for i=length(value):-1:1
        C(i) = {value(i)};
    end
end