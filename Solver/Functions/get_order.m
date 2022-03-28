function order = get_order(value)
    % Get order of magnitude of a number in base 10
    %
    % Args:
    %     value (float): number in base 10
    %
    % Returns:
    %     order (float): order of magnitude of a number in base 10
    
    order = floor(log(abs(value)) ./ log(10));
end