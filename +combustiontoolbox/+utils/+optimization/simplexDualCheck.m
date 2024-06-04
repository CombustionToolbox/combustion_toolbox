function [x, x_min] = simplexDualCheck(A_eq, b_eq, c, A_ineq, b_ineq)
    % Use Matlab's simplex method to solve the next linear programming problem:
    %     * max(min x) -> max t,
    %     * A * x = b,
    %     *    x >= 0.
    % 
    % Args:
    %    A_eq (float): Coefficient matrix for the equality constraints (A * x = b)
    %    b_eq (float): Right-hand side vector of the equality constraints
    %    c (float): Coefficient vector of the objective function max t -> min -t
    %    A_ineq (float): Coefficient matrix for the inequality constraints (x_j - t >= 0)
    %    b_ineq (float): Right-hand side vector of the inequality constraints
    %
    % Returns:
    %    x (float): Optimal solution vector
    %    x_min (float): Minimum value of x
    %
    % Example:
    %    [x, x_min] = simplexDualCheck(A_eq, b_eq, c, A_ineq, b_ineq)

    % Definitions
    [~, n] = size(A_ineq);
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');

    % Solve the optimization problem
    [x, ~, ~, ~] = linprog(c, A_ineq, b_ineq, A_eq, b_eq, zeros(n, 1), [], options);

    % Remove t
    x_min = x(end);
    x(end) = [];
end