function [listSpecies, moles] = getAir(FLAG_IDEAL_AIR)
    % Define the chemical species that compound air and its molar composition
    %
    % Args:
    %     FLAG_IDEAL_AIR (bool): Flag indicating consider ideal or non-ideal air mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     * listSpecies (cell): List of species
    %     * moles (float): Array with molar composition
    
    % Ideal air
    if FLAG_IDEAL_AIR
        listSpecies = {'N2', 'O2'};
        moles = [79, 21] / 21;
        return
    end
    
    % Air
    listSpecies = {'N2', 'O2', 'Ar', 'CO2'};
    moles = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
end
