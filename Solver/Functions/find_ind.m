function ind = find_ind(S,NameSpecies)
for i=length(S):-1:1
    ind(i) = find(strcmp(NameSpecies,S{i}));
end
