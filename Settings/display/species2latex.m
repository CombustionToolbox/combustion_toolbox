function speciesLatex = species2latex(species)
    % Convert string of a species into latex

    % Check if is a condensed species
    index = strfind(species, 'bgrb');
    if ~isempty(index)
        species = strcat(species(1:index-1), '$_{(gr)}$', species(index+4:end));
    else
        index = strfind(species, 'bLb');
        if ~isempty(index)
            species = strcat(species(1:index-1), '$_{(L)}$', species(index+3:end));
        end
    end
    % Check numbers
    index = regexp(species, '[0-9]');
    N = length(index);
    if ~N
        speciesLatex = species;
        return;
    end
    speciesLatex = [];
    pos1 = 1;
    for i=1:N
        pos2 = index(i) - 1;
        speciesLatex = strcat(speciesLatex, species(pos1:pos2), '$_', species(pos2 + 1), '$');
        pos1 = pos2 + 2;
    end
    if pos1 < length(species)
        if isletter(species(pos1))
            aux = 0;
        else
            aux = 1;
        end
        speciesLatex = [speciesLatex, species(pos1+aux:end)];
    end
end