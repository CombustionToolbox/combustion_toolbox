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

    for i = length(listSpecies):-1:1
        species = listSpecies{i};

        try
            value(1, i) = DB.(species).(property);
        catch
            value(1, i) = 0;
        end

    end

end
