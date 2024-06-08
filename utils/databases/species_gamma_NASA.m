function gamma = species_gamma_NASA(species, T, DB)
    % Compute adiabatic index of the species [-] at the given temperature
    % [K] using piecewise cubic Hermite interpolating polynomials and
    % linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Range of temperatures to evaluate [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     gamma (float): Adiabatic index [-]
    %
    % Example:
    %     gamma = species_gamma_NASA('H2O', 300:100:6000, DB)

    [cp, cv] = species_cP_NASA(species, T, DB);
    gamma = cp / cv;

    assert(~isnan(gamma), 'Adibatic index equal NaN');
end