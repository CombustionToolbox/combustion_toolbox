function cV = species_cV(species, T, DB)
    % Compute specific heat at constant volume [J/(mol-K)] of the species
    % at the given temperature [K] using piecewise cubic Hermite
    % interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cV (float): Specific heat at constant volume [J/(mol-K)]

    try
        cV = DB.(species).cPcurve(T) - 8.31446261815324;
    catch
        cV = 0; % cryogenic species
    end

end
