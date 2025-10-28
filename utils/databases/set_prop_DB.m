function value = set_prop_DB(listSpecies, property, DB)
    % Function that gets the vector of the defined property for the given
    % set of species
    %
    % Args:
    %     listSpecies (cell): List of species
    %     property (str): Property to obtain from the database
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     value (float): Property vector
    %
    % Example:
    %     value = set_prop_DB({'H2O', 'CO2'}, 'hf', DB)

    % Preallocate for better performance
    numSpecies = length(listSpecies);
    value = zeros(1, numSpecies);
    
    for i = 1:numSpecies
        species = listSpecies{i};

        try
            value(i) = DB.(species).(property);
        catch
            value(i) = 0;
        end

    end

end
