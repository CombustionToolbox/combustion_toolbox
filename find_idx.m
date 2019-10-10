function idx = find_idx(S,NameSpecies)
for i=length(S):-1:1
    idx(i) = find(strcmp(NameSpecies,S{i}));
end
