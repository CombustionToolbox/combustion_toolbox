function x = simplex(A, b, c)
    % Use revised simplex method to solve the linear programming problem:
    %     * min c' * x,
    %     * A * x = b,
    %     *    x >= 0
    %
    % Args:
    %    A (float): Coefficient matrix for the equality constraints (A * x = b)
    %    b (float): Right-hand side vector of the equality constraints
    %    c (float): Coefficient vector of the objective function (minimize c' * x)
    %
    % Returns:
    %    x (float): Optimal nonnegative solution vector
    %
    % Example:
    %    x = simplex(A, b, c)

    % Tolerance
    tol = 1e-4;

    % Normalize shapes
    b = b(:);
    c = c(:);
    [m, n] = size(A);

    assert(length(b) == m, 'simplex:DimensionMismatch', 'b length must match rows of A.');
    assert(length(c) == n, 'simplex:DimensionMismatch', 'c length must match columns of A.');

    itMax = max(100, 20 * (m + n));

    % Remove zero rows and normalize rows so b is nonnegative
    rowNorm = max(abs(A), [], 2);
    FLAG_ZERO_ROW = rowNorm <= tol;

    if any(abs(b(FLAG_ZERO_ROW)) > tol)
        error('simplex:Infeasible', 'Zero constraint row has nonzero right-hand side.');
    end

    A(FLAG_ZERO_ROW, :) = [];
    b(FLAG_ZERO_ROW) = [];
    rowNorm(FLAG_ZERO_ROW) = [];

    FLAG_NEGATIVE_B = b < 0;
    A(FLAG_NEGATIVE_B, :) = -A(FLAG_NEGATIVE_B, :);
    b(FLAG_NEGATIVE_B) = -b(FLAG_NEGATIVE_B);

    rowScale = max(rowNorm, 1);
    A = A ./ rowScale;
    b = b ./ rowScale;

    % Remove linearly dependent rows
    [A, b] = independentRows(A, b, tol);
    [m, n] = size(A);

    if m == 0
        x = zeros(n, 1);
        return
    end

    % Phase I: add artificial variables to obtain an explicit feasible basis
    A_phase = [A, eye(m)];
    c_phase = [zeros(n, 1); ones(m, 1)];
    basis = n + (1:m);

    [x_phase, basis] = solvePhase(A_phase, b, c_phase, basis, tol, itMax);
    phase1Objective = c_phase' * x_phase;

    if phase1Objective > max(1, norm(b, inf)) * tol
        error('simplex:Infeasible', 'Phase I could not find a feasible solution.');
    end

    basis = removeArtificialBasis(A_phase, A, basis, n, tol);

    % Phase II: optimize original objective from the feasible basis
    x = solvePhase(A, b, c, basis, tol, itMax);

    x(abs(x) <= tol) = 0;
end

function [A, b] = independentRows(A, b, tol)
    % Keep a full-row-rank subset selected by QR pivoting

    % Definitions
    [m, ~] = size(A);

    % Check empty constraint matrix
    if m == 0
        return
    end

    % Estimate row rank from column pivoting of A'
    [~, R, rowOrder] = qr(A', 0);
    diagR = abs(diag(R));

    % Get numerical rank
    if isempty(diagR)
        rankA = 0;
    else
        rankA = sum(diagR > max(size(A)) * eps(max(diagR)));
    end

    % Return if rows are already independent
    if rankA == m
        return
    end

    % Check consistency of rank-zero constraints
    if rankA == 0
        if norm(b, inf) > tol
            error('simplex:Infeasible', 'Rank-deficient constraints are inconsistent.');
        end

        A = zeros(0, size(A, 2));
        b = zeros(0, 1);
        return
    end

    % Keep independent rows selected by QR pivoting
    indexKeep = sort(rowOrder(1:rankA));
    A_keep = A(indexKeep, :);
    b_keep = b(indexKeep);

    % Check consistency before discarding dependent rows
    if norm(A * pinv(A_keep) * b_keep - b, inf) > max(1, norm(b, inf)) * tol
        error('simplex:Infeasible', 'Rank-deficient constraints are inconsistent.');
    end

    % Update reduced system
    A = A_keep;
    b = b_keep;
end

function basis = removeArtificialBasis(A_phase, A, basis, numOriginalVariables, tol)
    % Pivot artificial variables out of the basis after Phase I

    % Definitions
    numRows = length(basis);

    % Loop through basis rows
    for row = 1:numRows
        % Skip original variables already in the basis
        if basis(row) <= numOriginalVariables
            continue
        end

        % Build current Phase I basis
        B = A_phase(:, basis);

        % Get original nonbasic candidate columns
        FLAG_NONBASIC_ORIGINAL = true(1, numOriginalVariables);
        FLAG_NONBASIC_ORIGINAL(basis(basis <= numOriginalVariables)) = false;
        candidates = find(FLAG_NONBASIC_ORIGINAL);
        pivotColumn = [];

        % Search for an original column that can replace the artificial one
        for candidate = candidates
            direction = B \ A(:, candidate);

            if abs(direction(row)) > tol
                pivotColumn = candidate;
                break
            end
        end

        % Check that Phase I found a removable artificial basis variable
        if isempty(pivotColumn)
            error('simplex:DegenerateBasis', 'Could not remove artificial variable from basis.');
        end

        % Replace artificial column by original column
        basis(row) = pivotColumn;
    end
end

function [x, basis] = solvePhase(A, b, c, basis, tol, itMax)
    % Dense revised simplex with maintained basis inverse

    % Definitions
    [m, n] = size(A);
    basis = basis(:).';
    FLAG_BASIC = false(1, n);
    FLAG_BASIC(basis) = true;

    % Initialize basis inverse and basic solution
    Binv = basisInverse(A, basis);
    xBasis = Binv * b;

    % Initialization
    it = 0;
    refactorEvery = 50;

    % Loop
    while true
        % Update iteration number
        it = it + 1;

        % Check iteration limit
        if it > itMax
            error('simplex:IterationLimit', 'Simplex iteration limit exceeded.');
        end

        % Periodically refactorize the basis inverse to control roundoff
        if mod(it, refactorEvery) == 0
            Binv = basisInverse(A, basis);
            xBasis = Binv * b;
        end

        % Remove numerical noise from the basic solution
        xBasis(abs(xBasis) <= tol) = 0;

        % Check primal feasibility of the current basis
        if any(xBasis < -tol)
            error('simplex:InvalidBasis', 'Current basis is primal infeasible.');
        end

        % Compute simplex multipliers
        lambda = Binv.' * c(basis);

        % Compute reduced costs for nonbasic variables
        nonBasis = find(~FLAG_BASIC);
        reducedCosts = c(nonBasis) - A(:, nonBasis)' * lambda;

        % Check optimal solution
        if all(reducedCosts >= -tol)
            x = zeros(n, 1);
            x(basis) = max(xBasis, 0);
            return
        end

        % Find entering variable with Bland's rule
        enteringPosition = find(reducedCosts < -tol, 1, 'first');
        entering = nonBasis(enteringPosition);

        % Compute pivot direction
        direction = Binv * A(:, entering);
        FLAG_POSITIVE = direction > tol;

        % Check unbounded problem
        if ~any(FLAG_POSITIVE)
            error('simplex:Unbounded', 'Linear program is unbounded.');
        end

        % Find leaving row using the minimum ratio test
        ratios = inf(m, 1);
        ratios(FLAG_POSITIVE) = xBasis(FLAG_POSITIVE) ./ direction(FLAG_POSITIVE);
        minRatio = min(ratios);

        % Break ties with the smallest-index basic variable
        leavingCandidates = find(abs(ratios - minRatio) <= tol);
        [~, indexTie] = min(basis(leavingCandidates));
        leavingRow = leavingCandidates(indexTie);
        leaving = basis(leavingRow);
        pivot = direction(leavingRow);

        % Check numerical pivot
        if abs(pivot) <= tol
            error('simplex:NumericalPivot', 'Pivot is too small.');
        end

        % Update basic solution
        theta = xBasis(leavingRow) / pivot;
        xBasis = xBasis - direction * theta;
        xBasis(leavingRow) = theta;

        % Update basis inverse with a Gauss-Jordan pivot update
        pivotRow = Binv(leavingRow, :) / pivot;
        Binv = Binv - direction * pivotRow;
        Binv(leavingRow, :) = pivotRow;

        % Update basis bookkeeping
        basis(leavingRow) = entering;
        FLAG_BASIC(leaving) = false;
        FLAG_BASIC(entering) = true;
    end
end

function Binv = basisInverse(A, basis)
    % Build the dense inverse of the current basis

    % Definitions
    m = length(basis);

    % Solve B * Binv = I
    Binv = A(:, basis) \ eye(m);
end
