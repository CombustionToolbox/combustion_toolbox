function index = find_ind(LS, species)
    % Find the index of the species based on the given list (LS)
    %
    % Args:
    %     LS (cell):      List of species
    %     species (cell): Species to find index values
    %
    % Returns:
    %     index (float):  List with the index of the species based on the given list (LS)

    if length(species) > 1
        NS = length(LS);
        if ischar(species)
            Nspecies = 1;
            species = {species};
        else
            Nspecies = length(species);
        end
        index = [];
        for j=1:Nspecies
            i = 0;
            while i < NS
                i = i + 1;
                if strcmp(species{j}, LS{i})
                    index = [index, i];
                    break
                end 
            end
        end
    else
        index = ismember(LS, species);
        index = find(index);
    end
end
