function index = get_index_ions(species)
    index = (contains(species, 'minus') | contains(species, 'plus')) & ~contains(species, 'cyclominus');
end