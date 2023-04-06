function cP = species_cP(species, T, DB)
    % Compute specific heat at constant pressure [J/(mol-K)] of the species
    % at the given temperature [K] using piecewise cubic Hermite
    % interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cP (float): Specific heat at constant pressure in molar basis [J/(mol-K)]
    %
    % Example:
    %     cP = species_cP('H2O', 300, DB)

    try
        cP = DB.(species).cPcurve(T);
    catch
        cP = 0; % cryogenic species
    end

end
