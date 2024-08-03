function x = simplex(A, b, c)
    % Use simplex method to solve the next linear programming problem:
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
    %    x = simplex(A, b, c)

    % Get initial feasible solution
    [A, b, m, n] = initialization(A, b);
    
    % Solve the optimization problem
    tab = simplexMethod(A, b, c);
    
    % Get the solution vector from the final tableau
    x = getSolution(tab, m, n);
end

function [A, b, m, n] = initialization(A, b)
    % Get initial feasible solution

    % Definitions
    [m, n] = size(A);
    I = eye(m);
    A0 = [A I];
    c0 = [-sum(A0(:, 1:n)), zeros(1, m)];

    % Get initial feasible solution
    tab = simplexMethod(A0, b, c0');
    
    % Update A and b to feasible solution
    A = tab(1:m, 1:n);
    b = tab(1:m, end);
end

function tab = simplexMethod(A, b, c)
    % Solve the linear programming problem using the Simplex method

    % Definitions
    [m, n] = size(A);
    tab = [A b; c' 0];
    tol = 1e-4;
    
    % Initialization
    it = 0;
    
    % Loop
    while true
        % Update iteration number
        it = it + 1;
        
        % Find the pivot column (most negative value in the last row)
        [min_val, pivot_col] = min(tab(end, 1:n)); 
        
        % Check optimal solution
        if min_val >= -tol
            break
        end

        % Find the pivot row using the minimum ratio test
        pivot_col_vals = tab(1:m, pivot_col);
        ratios = tab(1:m, end) ./ pivot_col_vals;
        ratios(pivot_col_vals <= tol) = Inf;
        [min_ratio, pivot_row] = min(ratios);

        % Check unbounded problem
        if isinf(min_ratio)
            error('The problem is unbounded.');
        end

        % Update the tableau using the pivot element
        pivot_element = tab(pivot_row, pivot_col);
        tab(pivot_row, :) = tab(pivot_row, :) / pivot_element;

        % Update the rows
        if n < 100
            rows_to_update = true(1, m + 1);
            rows_to_update(pivot_row) = false;
            tab(rows_to_update, :) = tab(rows_to_update, :) - tab(rows_to_update, pivot_col) * tab(pivot_row, :);
            continue
        end

        for i = [1:pivot_row-1, pivot_row+1:m+1]
            tab(i, :) = tab(i, :) - tab(i, pivot_col) * tab(pivot_row, :);
        end

    end

end

function x = getSolution(tab, m, n)
    % Get the solution vector from the final tableau
    
    % Definitions
    FLAG = tab(1:m, 1:n) ~= 0;
    tol = 1e-4;
    
    % Initialization
    x = zeros(n, 1);

    % Get the indices of the basic variables
    indexBasic = find(sum(FLAG) == 1 & sum(tab(1:m, 1:n) == 1) & tab(end, 1:n) < tol);

    % Assign values to the basic variables in the solution vector
    x(indexBasic) = FLAG(:, indexBasic)' * tab(1:m, end);
end