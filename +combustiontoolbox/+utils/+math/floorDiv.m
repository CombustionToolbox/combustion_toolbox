function x = floorDiv(a, b)
    % Round the result of division toward negative infinity
    %
    % Args:
    %     a (float): Value a
    %     b (float): Value b
    %
    % Returns:
    %     x (float): Floor division
    %
    % Example:
    %     x = floorDiv(5, 2);

    x = floor(a ./ b);
end