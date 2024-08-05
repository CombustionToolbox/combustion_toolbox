function g0 = species_g0_NASA(obj, species, temperature)
    % Compute Compute Gibbs energy [J/mol] of the species at the given
    % temperature [K] using NASA's 9 polynomials
    %
    % Args:
    %     obj (NasaDatabase): NasaDatabase object
    %     species (char): Chemical species
    %     temperature (float): Range of temperatures to evaluate [K]
    %
    % Returns:
    %     g0 (float): Gibbs energy in molar basis [J/mol]
    %
    % Example:
    %     g0 = species_g0_NASA(NasaDatabase(), 'H2O', 300:100:6000)

    % Definitions
    R0 = combustiontoolbox.common.Constants.R0; % Universal Gas Constant [J/(mol-K)];

    % Unpack NASA's polynomials coefficients
    [a, b, ~, tExponents, ctTInt] = obj.getCoefficients(species, obj.species);
    
    % Compute specific enthalpy [J/mol]
    for i = length(temperature):-1:1
        T = temperature(i);

        if ctTInt > 0
            % Compute interval temperature
            tInterval = obj.getIndexTempereratureInterval(species, T, obj.species);

            % Compute Gibbs energy
            g0(i) = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1/2, 1 + log(T), 1 - log(T), -1/2, -1/6, -1/12, -1/20, 0]) + b{tInterval}(1) / T - b{tInterval}(2));
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
        else
            g0(i) = obj.species.(species).Hf0;
        end

    end

end