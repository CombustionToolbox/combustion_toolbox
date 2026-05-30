function [x, x_min] = simplexDual(A, b)
    % Use simplex method to solve the next linear programming problem:
    %     * max(min x) -> max t,
    %     * A * x = b,
    %     *    x >= 0.
    % 
    % Args:
    %    A (float): Coefficient matrix for the equality constraints (A * x = b)
    %    b (float): Right-hand side vector of the equality constraints
    %
    % Returns:
    %    Tuple containing
    %
    %    * x (float): Optimal solution vector
    %    * x_min (float): Minimum value of x
    %
    % Example:
    %    [x, x_min] = simplexDual(A, b)

    % Definitions
    [~, n] = size(A);
    b = b(:);
    sumA = sum(A, 2);

    % Coefficients for the objective function max t -> min -t
    c = zeros(n + 1, 1);
    c(end) = -1;

    % Set equality constraints (A * z + sum(A, 2) * t = b)
    A_eq = [A, sumA];

    % Solve the linear programming problem using the simplex method
    zt = combustiontoolbox.utils.optimization.simplex(A_eq, b, c);

    % Recover minimum value
    x_min = max(zt(end), 0);

    % Recover solution (x = z + t)
    x = zt(1:n) + x_min;
end
