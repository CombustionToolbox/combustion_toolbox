function dxdy = compute_first_derivative(x, y)
    % Compute first central derivate using a non-uniform grid
    %
    % Args:
    %     x (float): Grid values
    %     y (float): Values for the corresponding grid
    %
    % Returns:
    %     dxdy (float): Value of the first derivate for the given grid and its corresponding values

    h = y(2:end) - y(1:end-1);
    hmax = max(h);
    mu = h / hmax;

    dxdy = zeros(1, length(h));

    dxdy(1) = ((x(2) - x(1)) ./ h(1));
    for i = 2:length(mu)-1
        dxdy(i) = (mu(i)^2 * x(i+1) - (mu(i)^2 - mu(i+1)^2) * x(i) - mu(i+1)^2 * x(i-1)) / ((mu(i)^2 * mu(i+1) + mu(i) * mu(i+1)^2) * hmax);
    end
    dxdy(end) = ((x(end) - x(end-1)) ./ h(end));
end