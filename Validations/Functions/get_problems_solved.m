function value = get_problems_solved(mix, variable)
    % Get problems solved based on the length of the given variable
    for i=length(mix):-1:1
        value(i) = length(mix{i}.(variable));
    end
    value = sum(value);
end