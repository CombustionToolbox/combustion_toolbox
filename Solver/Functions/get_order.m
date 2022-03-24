function order = get_order(value)
    % Get order of magnitude of a number in base 10
    order = floor(log(abs(value)) ./ log(10));
end