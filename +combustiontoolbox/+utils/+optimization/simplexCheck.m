function x = simplexCheck(A, b, c)
    % Use Matlab's simplex method to solve the next linear programming problem:
    %     * min c'x,
    %     * A * x = b,
    %     *    x >= 0.
    % 
    % Args:
    %    A (float): Coefficient matrix for the equality constraints (Ax = b)
    %    b (float): Right-hand side vector of the equality constraints
    %    c (float): Coefficient vector of the objective function (minimize c'x)
    %
    % Return:
    %    x (float): Optimal solution vector
    %
    % Example:
    %    x = simplexCheck(A, b, c)

    % Definitions
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'off');
    
    % Solve the optimization problem
    [x, ~, ~, ~] = linprog(c, [], [], A, b', zeros(length(c), 1), [], options);
end