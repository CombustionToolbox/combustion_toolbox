function [cP, cV] = species_cP_NASA(species, temperature, DB)
    % Compute specific heats at constant pressure and at constant volume
    % [J/(mol-K)] of the species at the given temperature [K] using NASA's
    % 9 polynomials
    %
    % Args:
    %     species (str): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     - cP (float): Specific heat at constant pressure [J/(mol-K)]
    %     - cV (float): Specific heat at constant volume   [J/(mol-K)]

    % Definitions
    R0 = 8.31446261815324; % Universal Gas Constant [J/(mol-K)];
    % Unpack NASA's polynomials coefficients
    [a, ~, tRange, tExponents, ctTInt] = unpack_NASA_coefficients(species, DB);
    % Compute specific heat at constant pressure
    for i = length(temperature):-1:1
        T = temperature(i);

        if DB.(species).ctTInt > 0
            % Compute interval temperature
            tInterval = compute_interval_NASA(species, T, DB, tRange, ctTInt);
            % Compute specific specific heat at constant pressure
            cP(i) = R0 * (sum(a{tInterval} .* T.^tExponents{tInterval}));
            cV(i) = cP(i) - R0;
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            cP(i) = 0;
            cV(i) = 0;
        end

    end

end
