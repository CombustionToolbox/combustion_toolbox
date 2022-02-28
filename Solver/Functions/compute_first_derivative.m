function dxdy = compute_first_derivative(x, y)
    % Compute first central derivate with non-uniform grid
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