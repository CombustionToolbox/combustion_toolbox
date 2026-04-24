function speciesArray = addPengRobinsonProperties(speciesArray)
    % Assign parameters for pure species for Peng-Robinson EoS.
    %
    % References:
    % [1] THERMODYNAMIC AND THERMOPHYSICAL PROPERTIES OF HUMID AIR BY USING CUBIC PENG-ROBINSON EOS
    % [2] PROPERTIES OF GASES, Isidoro Martinez

    % 1. Initialize all existing species in the struct with NaN.
    % This prevents [] (empty) values which would crash vectorized logical checks.
    fnames = fieldnames(speciesArray);
    for i = 1:length(fnames)
        speciesArray.(fnames{i}).Tcritical = NaN;
        speciesArray.(fnames{i}).Pcritical = NaN;
        speciesArray.(fnames{i}).acentricFactor = NaN;
    end

    % 2. Assign known values
    speciesArray = addSpecies(speciesArray, 'Ar', 150.87, 48.98, 0.000); 
    speciesArray = addSpecies(speciesArray, 'CH4', 190.56, 45.99, 0.011);
    speciesArray = addSpecies(speciesArray, 'C2H4', 282.35, 50.42, 0.087);  
    speciesArray = addSpecies(speciesArray, 'C2H6', 305.32, 48.72, 0.099);
    speciesArray = addSpecies(speciesArray, 'C3H8', 369.83, 42.48, 0.152);
    speciesArray = addSpecies(speciesArray, 'CO', 132.86, 34.94, 0.045); 
    speciesArray = addSpecies(speciesArray, 'CO2', 304.13, 73.75, 0.225);
    speciesArray = addSpecies(speciesArray, 'O2', 154.581, 50.429, 0.0222);  
    speciesArray = addSpecies(speciesArray, 'N2', 126.19, 33.96, 0.037);  
    speciesArray = addSpecies(speciesArray, 'NO', 180.15, 64.80, 0.588);
    speciesArray = addSpecies(speciesArray, 'NO2', 431.20, 101.32, 0.851);
    speciesArray = addSpecies(speciesArray, 'H2', 33.19,  13.13, -0.219);
    speciesArray = addSpecies(speciesArray, 'H2O', 647.14, 220.64, 0.344);  
    speciesArray = addSpecies(speciesArray, 'H2ObLb', 647.14, 220.64, 0.344);
    speciesArray = addSpecies(speciesArray, 'C12H10', 789.26, 38.47, 0.3720);
end

% SUB-PASS FUNCTIONS
function speciesArray = addSpecies(speciesArray, species, Tcritical, Pcritical, acentricFactor)
    % Only assign if the species name actually exists as a field in the struct
    assert(isfield(speciesArray, species), 'Species "%s" not found in the input struct. Check for typos or missing entries.', species);

    speciesArray.(species).Tcritical = Tcritical;           % [K]
    speciesArray.(species).Pcritical = Pcritical;           % [bar]
    speciesArray.(species).acentricFactor = acentricFactor; % [-]
end