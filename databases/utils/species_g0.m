function g0 = species_g0(species, T, DB)
    % Compute Gibbs energy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     g0 (float): Gibbs energy [kJ/mol]

    try
        g0 = DB.(species).g0curve(T) / 1000;
    catch
        g0 = DB.(species).g0 / 1000; % cryogenic species
    end

end
