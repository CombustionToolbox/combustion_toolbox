function gamma = species_gamma(species, T, DB)
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

    gamma = species_cP(species, T, DB) / species_cV(species, T, DB);

    assert(~isnan(gamma), 'Adibatic index equal NaN');
end