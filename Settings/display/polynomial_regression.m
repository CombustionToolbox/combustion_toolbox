function y_poly = polynomial_regression(x, y, n)
    % Obtain polynomial regression for the given dataset (x, y)
    % and polynomial order
    %
    % Args:
    %     x (float): x values
    %     y (float): y values
    %     n (float): polynomial order
    % Returns:
    %     y_poly (float): y values of the polynomial regression

    % Obtain polynomial fit
    p = polyfit(x, y, n);
    % Evaluate polynomial
    y_poly = polyval(p, x);
end
