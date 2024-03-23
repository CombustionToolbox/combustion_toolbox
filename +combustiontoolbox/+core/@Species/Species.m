classdef Species

    properties
        name          % Name chemical species
        fullname      % Fullname chemical species
        comments      % Additional comments from database
        formula       % Chemical formula
        W             % Weight mass
        hf            % Enthalpy of formation
        ef            % Internal energy of formation
        phase         % Phase
        T             % Temperature
        Tintervals    % Number temperature intervals 
        Trange        % Polynomials temperature intervals
        Texponents    % Exponents polynomials
        a             % Coefficients a polynomials
        b             % Coefficients b polynomials
        cpcurve       % Gridded interpolant object with specific heat at constant pressure of the individual species
        h0curve       % Gridded interpolant object with enthalpy of the individual species
        s0curve       % Gridded interpolanb object with entropy of the individual species
        g0curve       % Gridded interpolanb object with Gibbs free energy of the individual species
        elementMatrix % Element matrix
    end

end