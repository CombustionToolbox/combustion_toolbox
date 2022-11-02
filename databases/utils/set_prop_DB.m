function value = set_prop_DB(LS, property, DB)
    % Function that gets the vector of the defined property for the given
    % set of species
    %
    % Args:
    %     LS (cell): List of species
    %     property (str): Property to obtain from the database
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     value (float): Property vector

    for i = length(LS):-1:1
        species = LS{i};

        try
            value(1, i) = DB.(species).(property);
        catch
            value(1, i) = 0;
        end

    end

end
