function g0 = species_g0_NASA(species, temperature, DB)
    % Compute Compute Gibbs energy [kJ/mol] of the species at the given
    % temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (str): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     g0 (float): Gibbs energy [kJ/mol]

    % Definitions
    R0 = 8.31446261815324; % Universal Gas Constant [J/(mol-K)];
    % Unpack NASA's polynomials coefficients
    [a, b, tRange, tExponents, ctTInt] = unpack_NASA_coefficients(species, DB);
    % Compute specific enthalpy [kJ/mol]
    for i = length(temperature):-1:1
        T = temperature(i);

        if DB.(species).ctTInt > 0
            % Compute interval temperature
            tInterval = compute_interval_NASA(species, T, DB, tRange, ctTInt);
            % Compute Gibbs energy
            % h0 = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1, log(T), 1, 1/2, 1/3, 1/4, 1/5, 0]) + b{tInterval}(1)/T);
            % s0 = R0 * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1/2, -1, log(T), 1, 1/2, 1/3, 1/4, 0]) + b{tInterval}(2));
            % g0(i) = h0 - T * s0;

            g0(i) = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1/2, 1 + log(T), 1 - log(T), -1/2, -1/6, -1/12, -1/20, 0]) + b{tInterval}(1) / T - b{tInterval}(2));
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            g0(i) = DB.(species).Hf0;
        end

    end

    % Change units [kJ/mol]
    g0 = g0 * 1e-3;
end
