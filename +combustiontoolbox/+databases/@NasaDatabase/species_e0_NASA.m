function [e0, DeT] = species_e0_NASA(species, temperature, DB)
    % Compute internal energy and the thermal internal energy [J/mol] of
    % the species at the given temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     species (char): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     * e0 (float): Internal energy in molar basis [J/mol]
    %     * DeT (float): Thermal internal energy in molar basis [J/mol]
    %
    % Example:
    %     [e0, DeT] = species_e0_NASA('H2O', 300:100:6000, DB)

    % Definitions
    R0 = combustiontoolbox.common.Constants.R0; % Universal Gas Constant [J/(mol-K)];
    hf0 = DB.(species).hf; % [J/mol];
    Tref = 298.15; % [K]

    % Unpack NASA's polynomials coefficients
    [a, b, tRange, tExponents, ctTInt, txFormula, swtCondensed] = unpack_NASA_coefficients(species, DB);

    % Get elements
    elements = combustiontoolbox.core.elements().setElements();

    % Get element matrix of the species
    elementMatrix = setElementMatrix(txFormula, elements);

    % Compute change in moles of gases during the formation reaction of a
    % mole of that species starting from the elements in their reference state
    Delta_n = compute_change_moles_gas_reaction(elementMatrix, swtCondensed);
    % Compute specific enthalpy [J/mol]
    for i = length(temperature):-1:1
        T = temperature(i);

        if ctTInt > 0
            % Compute interval temperature
            tInterval = obj.getIndexTempereratureInterval(species, T, obj.species);

            % Compute specific enthalpy from NASA's 9 polynomials
            h0(i) = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1 log(T) 1 1/2 1/3 1/4 1/5 0]) + b{tInterval}(1) / T);
            ef0 = hf0 - Delta_n * R0 * Tref;
            e0(i) = (ef0 + (h0(i) - hf0) - (1 - swtCondensed) * R0 * (T - Tref));
            DeT(i) = e0(i) - ef0;
            
        % If the species is only a reactant determine it's reference temperature
        % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
        % K. For cryogenic liquids, assigned enthalpies are given at their boiling
        % points instead of 298.15 K
        else
            Tref = tRange(1);
            e0(i) = hf0 - Delta_n * R0 * Tref;
            DeT(i) = 0;
        end

    end

    % Change units [kJ/mol]
    % e0 = e0 * 1e-3;
    % DeT = DeT * 1e-3;
end
