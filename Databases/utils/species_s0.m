function s0 = species_s0(species, T, DB)
    % Compute entropy [kJ/(mol-K)] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (str): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     s0 (float): Entropy [kJ/(mol-K)]

    try
        s0 = DB.(species).s0curve(T) / 1000;
    catch
        s0 = 0; % cryogenic species
    end

end
