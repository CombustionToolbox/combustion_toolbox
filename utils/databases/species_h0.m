function h0 = species_h0(species, T, DB)
    % Compute enthalpy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     h0 (float): Enthalpy in molar basis [kJ/mol]
    %
    % Example:
    %     h0 = species_h0('H2O', 300, DB)

    try
        h0 = DB.(species).h0curve(T) / 1000;
    catch
        h0 = DB.(species).h0 / 1000; % cryogenic species
    end

end
