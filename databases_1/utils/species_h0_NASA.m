function [h0, DhT] = species_h0_NASA(species, temperature, DB)
    % Compute enthalpy and thermal enthalpy [kJ/mol] of the species at the
    % given temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (str): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     - h0 (float): Enthalpy [kJ/mol]
    %     - DhT (float): Thermal enthalpy [kJ/mol]

    % Definitions
    R0 = 8.31446261815324; % Universal Gas Constant [J/(mol-K)];
    hf0 = DB.(species).hf; % [J/mol];
    % Unpack NASA's polynomials coefficients
    [a, b, tRange, tExponents, ctTInt] = unpack_NASA_coefficients(species, DB);
    % Compute specific enthalpy [kJ/mol]
    for i = length(temperature):-1:1
        T = temperature(i);

        if DB.(species).ctTInt > 0
            % Compute interval temperature
            tInterval = compute_interval_NASA(species, T, DB, tRange, ctTInt);
            % Compute specific enthalpy from NASA's 9 polynomials
            h0(i) = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1 log(T) 1 1/2 1/3 1/4 1/5 0]) + b{tInterval}(1) / T);
            DhT(i) = h0(i) - hf0;
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            h0(i) = hf0;
            DhT(i) = 0;
        end

    end

    % Change units [kJ/mol]
    h0 = h0 * 1e-3;
    DhT = DhT * 1e-3;
end
