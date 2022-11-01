function name = FullName2name(species)
    % Get full name of the given species
    %
    % Args:
    %     species (str): Chemical species
    %
    % Return:
    %     name (str): Full name of the given species

    FLAG_MILLENIUM = false;

    if contains(species, '_M')
        species = strrep(species, '_M', '');
        FLAG_MILLENIUM = true;
    end

    name = species;

    if isempty(name)
        return
    end

    if name(end) == '+'
        name = [name(1:end - 1) 'plus'];
    elseif name(end) == '-'
        name = [name(1:end - 1) 'minus'];
    end

    ind = regexp(name, '[()]');
    name(ind) = 'b';
    ind = regexp(name, '[.,+-]');
    name(ind) = '_';

    if regexp(name(1), '[0-9]')
        name = ['num_' name];
    end

    ind = regexp(name, '\x27');
    name(ind) = '_';

    if FLAG_MILLENIUM
        name = strcat(name, '_M');
    end

end
