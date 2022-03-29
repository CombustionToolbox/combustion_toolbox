function index = get_index_ions(species)
    % Get index of ions for the given list of species
    %
    % Args:
    %     species (str): List of species
    % 
    % Returns:
    %     index (float): Index of ions
    
    index = (contains(species, 'minus') | contains(species, 'plus')) & ~contains(species, 'cyclominus');
end