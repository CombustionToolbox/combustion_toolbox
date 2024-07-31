function [N, NP] = equilibriumGuess(N, NP, A0, muRT, b0, index, indexGas, indexIons, NG, molesGuess)
    % Initialize molar composition from a previous calculation or using the simplex algorithm
    %
    % Args:
    %     N (float): Mixture composition [mol]
    %     NP (float): Total moles of gaseous species [mol]
    %     A0 (float): Stoichiometric matrix
    %     muRT (float): Dimensionless chemical potentials
    %     b0 (float): Moles of each element
    %     index (float): Index of chemical species
    %     indexGas (float): Index of gaseous species
    %     indexIons (float): Index of ionized species
    %     NG (float): Number of gaseous species
    %     molesGuess (float): Mixture composition [mol] of a previous computation
    %
    % Returns:
    %     Tuple containing
    %
    %     * N (float): Mixture composition [mol]
    %     * NP (float): Total moles of gaseous species [mol]
    %
    % Example:
    %     [N, NP] = equilibriumGuess(N, NP, A0, muRT, b0, index, indexGas, indexIons, NG, molesGuess)
    
    % Get molar composition from a previous calculation
    if ~isempty(molesGuess)
        N(indexGas) = molesGuess(indexGas);
        NP = sum(molesGuess(indexGas));
        return
    end
    
    try
        % Get molar composition using the simplex method
        N = getSimplex(N, A0, muRT, b0, index, indexIons, NG);
    catch
        % Get molar composition using a uniform distribution
        N(index) = NP/NG;
    end

    % Recompute mol gaseous species
    NP = sum(N(indexGas));
end

% SUB-PASS FUNCTIONS
function N = getSimplex(N, A0, muRT, b0, index, indexIons, NG)
    % Get molar composition using the simplex algorithm

    % Definitions
    alpha = 0.01;
    FLAG_MINOR = true;

    % Initialization
    Nmin = 1e-2;
    Nminor = 0 * N;

    % Get major species
    Nmajor = combustiontoolbox.utils.optimization.simplex(A0, b0', muRT);

    % Remove ionized species from Nmajor
    Nmajor(indexIons) = 0;

    % Get minor species
    if FLAG_MINOR
        FLAG_MAXMIN = Nmajor > 0;
        indexPass = unique([1:NG, NG + index(FLAG_MAXMIN(NG + 1:end))]);
        [Nminor(indexPass), Nmin] = combustiontoolbox.utils.optimization.simplexDual(A0(:, indexPass), b0');
    end
    
    % Merge solutions
    N(index) = (1 - alpha) * Nmajor +  alpha * Nminor(index);
    N(index(N(index) == 0)) = alpha * Nmin;

    % Check
    % NmajorCheck = combustiontoolbox.utils.optimization.simplexCheck(A0, b0, muRT);

    % Check equality constrain (A * x = b)
    % error = A0 * N(index) - b0';
end