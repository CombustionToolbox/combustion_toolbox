function [index, NS] = check_temperature_range(system, T, index, NS, FLAG)
    % Remove species indeces out of the temperature range if FLAG = true,
    % e.g., linear extrapolation of the polynomial fits is not allowed.
    %
    % Args:
    %     system (ChemicalSystem): 
    %     T (float): Temperature
    %     ind (float): Vector with the species indeces
    %     NS (float): Number of species
    %     FLAG (bool): Flag indicating linear extrapolation of the polynomials fits
    %
    % Returns:
    %     Tuple containing
    %
    %     * ind (float): Updated vector with the species indeces
    %     * NS (float): Update number of species
    
    % Check FLAG
    if FLAG
        return
    end
    
    % Remove species outside the bounds
    for i = NS:-1:1
        species = system.listSpecies{index(i)};
        if T < system.species.(species).T(1) || T > system.species.(species).T(end)
            index(i) = [];
            NS = NS - 1;
        end

    end

end