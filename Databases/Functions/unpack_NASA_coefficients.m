function [a, b, tRange, tExponents, ctTInt, txFormula, swtCondensed] = unpack_NASA_coefficients(species, DB)
    % Unpack NASA's polynomials coefficients from database
    %
    %   Args:
    %     species (str) : Species
    %     DB (struct) : Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     - a (cell):
    %     - b (cell): 
    %     - tRange (cell): 
    %     - tExponents (cell): 
    %     - ctTInt (cell): 
    %     - txFormula (str): 
    %     - swtCondensed (float): 


    a = DB.(species).a;
    b = DB.(species).b;
    tRange = DB.(species).tRange;
    tExponents = DB.(species).tExponents;
    ctTInt = DB.(species).ctTInt;
    txFormula = DB.(species).txFormula;
    swtCondensed = DB.(species).swtCondensed;
end