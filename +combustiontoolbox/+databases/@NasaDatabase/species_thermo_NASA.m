function [cp0, cv0, h0, DhT, e0, DeT, s0, g0] = species_thermo_NASA(species, temperature, DB)
    % Compute thermodynamic function using NASA's 9 polynomials
    %
    % Args:
    %     species (char): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     * cP  (float): Specific heat at constant pressure in molar basis [J/(mol-K)]
    %     * cV  (float): Specific heat at constant volume in molar basis   [J/(mol-K)]
    %     * h0  (float): Enthalpy in molar basis [J/mol]
    %     * DhT (float): Thermal enthalpy in molar basis [J/mol]
    %     * e0  (float): Internal energy in molar basis [J/mol]
    %     * DeT (float): Thermal internal energy in molar basis [J/mol]
    %     * s0  (float): Entropy in molar basis [J/(mol-K)]
    %     * g0  (float): Gibbs energy in molar basis [J/mol]
    %
    % Example:
    %     [cP, cV, h0, DhT, e0, DeT, s0, g0] = species_thermo_NASA('H2O', 300:100:6000, DB)

    % Definitions
    R0 = combustiontoolbox.common.Constants.R0; % Universal Gas Constant [J/(mol-K)];
    hf0 = DB.(species).hf; % [J/mol];
    Tref = 298.15; % [K]
    % Unpack NASA's polynomials coefficients
    [a, b, tRange, tExponents, ctTInt, txFormula, phase] = unpack_NASA_coefficients(species, DB);
    % Get elements
    elements = set_elements();
    % Get element matrix of the species
    element_matrix = set_element_matrix(txFormula, elements);
    % Compute change in moles of gases during the formation reaction of a
    % mole of that species starting from the elements in their reference state
    Delta_n = compute_change_moles_gas_reaction(element_matrix, phase);
    % Compute thermodynamic properties
    for i = length(temperature):-1:1
        T = temperature(i);

        if DB.(species).ctTInt > 0
            % Compute interval temperature
            tInterval = compute_interval_NASA(species, T, DB, tRange, ctTInt);
            % Compute thermodynamic function
            cp0(i) = R0 * (sum(a{tInterval} .* T.^tExponents{tInterval}));
            h0(i) = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1 log(T) 1 1/2 1/3 1/4 1/5 0]) + b{tInterval}(1) / T);
            s0(i) = R0 * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1/2 -1 log(T) 1 1/2 1/3 1/4 0]) + b{tInterval}(2));
            ef0 = hf0 - Delta_n * R0 * Tref;
            e0(i) = (ef0 + (h0(i) - hf0) - (1 - phase) * R0 * (T - Tref));
            cv0(i) = cp0(i) - R0;
            DhT(i) = h0(i) - hf0;
            DeT(i) = e0(i) - ef0;
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            Tref = tRange(1);
            cp0 = 0;
            cv0 = 0;
            h0(i) = hf0;
            e0(i) = hf0 - Delta_n * R0 * Tref;
            g0(i) = DB.(species).Hf0;
            s0(i) = 0;
            DhT(i) = 0;
            DeT(i) = 0;
        end

    end

    % Change units
    % h0 = h0 * 1e-3; % [kJ/mol]
    % DhT = DhT * 1e-3; % [kJ/mol]
    % e0 = e0 * 1e-3; % [kJ/mol]
    % DeT = DeT * 1e-3; % [kJ/mol]
    % g0 = g0 * 1e-3; % [kJ/mol]
    % s0 = s0 * 1e-3; % [kJ/(mol-K)]
end