function ind = find_ind(S, species)
    % Find index of the species in the list S
    NS = length(S);
    i = 0;
    while i < NS
        i = i + 1;
        if strcmp(species, S{i})
            ind = i;
            return
        end
        
    end
    
    ind = false;
end
