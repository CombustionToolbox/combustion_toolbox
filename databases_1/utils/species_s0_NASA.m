function s0 = species_s0_NASA(species, temperature, DB)
    % Compute entropy [kJ/(mol-K)] of the species at the given temperature [K]
    % using NASA's 9 polynomials
    %
    % Args:
    %     species (str): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     s0 (float): Entropy [kJ/(mol-K)]

    % Definitions
    R0 = 8.31446261815324; % Universal Gas Constant [J/(mol-K)];
    % Unpack NASA's polynomials coefficients
    [a, b, tRange, tExponents, ctTInt] = unpack_NASA_coefficients(species, DB);
    % Compute entropy
    for i = length(temperature):-1:1
        T = temperature(i);

        if DB.(species).ctTInt > 0
            % Compute interval temperature
            tInterval = compute_interval_NASA(species, T, DB, tRange, ctTInt);
            % Compute entropy
            s0(i) = R0 * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1/2 -1 log(T) 1 1/2 1/3 1/4 0]) + b{tInterval}(2));
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            s0(i) = 0;
        end

    end

    % Change units [kJ/(mol-K)]
    s0 = s0 * 1e-3;
end
