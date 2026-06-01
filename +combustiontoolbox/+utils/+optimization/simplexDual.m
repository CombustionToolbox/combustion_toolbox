function [x, x_min] = simplexDual(A, b)
    % Use simplex method to solve the next linear programming problem:
    %     * max(min x) -> max t,
    %     * A * x = b,
    %     *    x >= 0
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
    A_original = A;
    b = b(:);
    [~, n] = size(A_original);

    % Compact max-min formulation with x = z + t, z >= 0, t >= 0
    A_eq = [A_original, sum(A_original, 2)];
    c = [zeros(n, 1); -1];

    % Solve equality-form linear program
    zt = combustiontoolbox.utils.optimization.simplex(A_eq, b, c);

    % Recover solution
    x_min = max(zt(end), 0);
    x = zt(1:n) + x_min;
end
