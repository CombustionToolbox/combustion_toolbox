function DhT = species_DhT(species, T, DB)
    % Compute thermal enthalpy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy [kJ/mol]

    try
        DhT = DB.(species).DhTcurve(T) / 1000;
    catch
        DhT = 0; % cryogenic species
    end

end
