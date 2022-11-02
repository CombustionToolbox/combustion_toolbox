function [a, b, tRange, tExponents, ctTInt, txFormula, phase] = unpack_NASA_coefficients(species, DB)
    % Unpack NASA's polynomials coefficients from database
    %
    % Args:
    %     species (str) : Chemical species
    %     DB (struct) : Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     * a (cell): Temperature coefficients
    %     * b (cell): Integration constants
    %     * tRange (cell): Ranges of temperatures [K]
    %     * tExponents (cell): Exponent coefficients
    %     * ctTInt (float): Number of intervals of temperatures
    %     * txFormula (str): Chemical formula
    %     * phase (float): 0 or 1 indicating gas or condensed phase, respectively

    a = DB.(species).a;
    b = DB.(species).b;
    tRange = DB.(species).tRange;
    tExponents = DB.(species).tExponents;
    ctTInt = DB.(species).ctTInt;
    txFormula = DB.(species).txFormula;
    phase = DB.(species).phase;
end
