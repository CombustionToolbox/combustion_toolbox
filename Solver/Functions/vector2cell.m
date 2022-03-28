function c = vector2cell(value)
    % Create cell array from vector
    %
    % Args:
    %     value (any): Vector with data of any type
    %
    % Returns:
    %     c (cell): Cell with the values of the vector

    for i=length(value):-1:1
        c(i) = {value(i)};
    end
end