function gamma = species_gamma_NASA(species, T, DB)
    % Compute adiabatic index of the species [-] at the given temperature
    % [K] using piecewise cubic Hermite interpolating polynomials and
    % linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     gamma (float): Adiabatic index [-]

    [cP, cV] = species_cP_NASA(species, T, DB);
    gamma = cP / cV;

    assert(~isnan(gamma), 'Adibatic index equal NaN');
end