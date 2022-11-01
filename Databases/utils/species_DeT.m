function DeT = species_DeT(species, T, DB)
    % Compute thermal internal energy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DeT (float): Thermal internal energy [kJ/mol]

    try
        DeT = DB.(species).DeTcurve(T) / 1000;
    catch
        DeT = 0; % cryogenic species
    end

end
