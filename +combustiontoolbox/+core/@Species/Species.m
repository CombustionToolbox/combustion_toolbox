classdef Species < handle
    % The :mat:func:`Species` class is used to store the properties of a chemical species.
    %
    % See also: :mat:func:`Database`, :mat:func:`NasaDatabase`

    properties
        name          % Name chemical species
        fullname      % Fullname chemical species
        refCode       % Reference date code
        comments      % Additional comments from database
        formula       % Chemical formula
        W             % Molecular weight [kg/mol]
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
    end

    properties (Hidden)
        elementMatrix % Element matrix
        FLAG_REFERENCE = false % Flag indicating if this species is a reference element/species
    end

    methods (Access = public)

        cp = getHeatCapacityPressure(obj, T)
        cv = getHeatCapacityVolume(obj, T)
        e0 = getInternalEnergy(obj, T)
        h0 = getEnthalpy(obj, T)
        s0 = getEntropy(obj, T)
        g0 = getGibbsEnergy(obj, T)
        DeT = getThermalInternalEnergy(obj, T)
        DhT = getThermalEnthalpy(obj, T)
        gamma = getAdiabaticIndex(obj, T)

        function elementMatrix = getElementMatrix(obj, elements)
            % Compute element matrix of the given species formula
            %
            % Args:
            %     elements (cell): List of elements
            %
            % Returns:
            %     elementMatrix(float): Element matrix. The first row refer to the index and the second to the number of atoms of each element contained in the species.
            %
            % Example:
            %     For CO2
            %
            %     elementMatrix = [7, 9; 1, 2]
            %
            %     That is, the species contains 1 atom of element 7 (C) and
            %     2 atoms of element 9 (O)
        
            % Definitions
            N = 40;
            NE = 5;
            step = 8;
            
            % Fill element matrix
            for i = 5:-1:1
                end0 = N - step * (NE - i);
                start0 = end0 - 5;
        
                element_i = obj.formula(start0 - 2:end0 - 6);
        
                if strcmp(element_i, '  ')
                    continue
                end
        
                if strcmp(element_i(2), ' ')
                    element_i = element_i(1);
                end
    
                elementMatrix(1, i) = find(strcmpi(elements, element_i));
                elementMatrix(2, i) = sscanf(obj.formula(start0:end0), '%f');
            end
        
        end

    end

end