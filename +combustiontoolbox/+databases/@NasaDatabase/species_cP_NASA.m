function [cp, cv] = species_cP_NASA(obj, species, temperature)
    % Compute specific heats at constant pressure and at constant volume
    % [J/(mol-K)] of the species at the given temperature [K] using NASA's
    % 9 polynomials
    %
    % Args:
    %     obj (NasaDatabase): NasaDatabase object
    %     species (char): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %
    % Returns:
    %     Tuple containing
    %
    %     * cp (float): Specific heat at constant pressure in molar basis [J/(mol-K)]
    %     * cv (float): Specific heat at constant volume in molar basis   [J/(mol-K)]
    %
    % Example:
    %     [cp, cv] = species_cP_NASA('H2O', 300:100:6000, DB)

    % Definitions
    R0 = combustiontoolbox.common.Constants.R0; % Universal Gas Constant [J/(mol-K)];

    % Unpack NASA's polynomials coefficients
    [a, ~, ~, tExponents, ctTInt] = obj.getCoefficients(species, obj.species);
    
    % Compute specific heat at constant pressure
    for i = length(temperature):-1:1
        T = temperature(i);

        if ctTInt > 0
            % Compute interval temperature
            tInterval = obj.getIndexTempereratureInterval(species, T, obj.species);

            % Compute specific specific heat at constant pressure
            cp(i) = R0 * (sum(a{tInterval} .* T.^tExponents{tInterval}));
            cv(i) = cp(i) - R0;
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            cp(i) = 0;
            cv(i) = 0;
        end

    end

end
