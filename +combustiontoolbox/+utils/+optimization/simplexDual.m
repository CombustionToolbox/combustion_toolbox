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
    %    x (float): Optimal solution vector
    %    x_min (float): Minimum value of x
    %
    % Example:
    %    [x, x_min] = simplexDual(A, b)

    % Definitions
    [m, n] = size(A);
    
    % Coefficients for the objective function max t -> min -t
    c = [zeros(n, 1); -1];
    
    % Set inequality constraints (x_j - t >= 0)
    A_ineq = [-eye(n), ones(n, 1)];
    b_ineq = zeros(n, 1);
    
    % Set equality constraints (A * x = b)
    A_eq = [A, zeros(m, 1)];
    b_eq = b;
    
    % Solve the linear programming problem using the simplex method
    x = combustiontoolbox.utils.optimization.simplex([A_ineq; A_eq], [b_ineq; b_eq], c);

    % Remove t
    x_min = x(end);
    x(end) = [];

    % Check
    % [x, x_min] = combustiontoolbox.utils.optimization.simplexDualCheck(A_eq, b_eq, c, A_ineq, b_ineq);
end