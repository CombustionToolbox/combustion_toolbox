classdef Species

    properties
        name          % Name chemical species
        fullname      % Fullname chemical species
        refCode       % Reference date code
        comments      % Additional comments from database
        formula       % Chemical formula
        W             % Molecular weight [g/mol]
        hf            % Enthalpy of formation at Tref from its reference species in their standard state [J/mol]
        hftoh0        % Enthalpy of formation at Tref relative to molar enthalpy at 0 K for standard state [J/mol]
        ef            % Internal energy of formation [J/mol]
        phase         % Phase
        T             % Temperature [K]
        Tref          % Temperature reference [K]
        Trange        % Temperature range intervals
        Tintervals    % Number of temperature intervals
        Texponents    % Exponents polynomials
        a             % Coefficients a polynomials
        b             % Coefficients b polynomials
        cpcurve       % Gridded interpolant object with specific heat at constant pressure of the individual species
        h0curve       % Gridded interpolant object with enthalpy of the individual species
        s0curve       % Gridded interpolant object with entropy of the individual species
        g0curve       % Gridded interpolant object with Gibbs free energy of the individual species
        elementMatrix % Element matrix
    end

end