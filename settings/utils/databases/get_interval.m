function tInterval = get_interval(species, T, DB)
    % Get interval of the NASA's polynomials from the Database (DB) for the
    % given species and temperature [K].
    %
    % Args:
    %     species (str): Chemical species
    %     T (float):     Temperature [K]
    %     DB (struct):   Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     tInterval (float): Index of the interval of temperatures

    for i = 1:DB.(species).ctTInt

        if (T >= DB.(species).tRange{i}(1)) && (T <= DB.(species).tRange{i}(2))
            break
        end

    end

    tInterval = i;
end
