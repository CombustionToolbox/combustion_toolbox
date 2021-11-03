function ind = find_ind(S, species)
    % Find index of the species in the list S
    for i=length(S):-1:1
        ind(i) = find(strcmp(species, S{i}));
    end
end
