names = fieldnames(strThProp);
%names = fieldnames(strMaster);
n_species = length(names);
for i = 1:n_species
    name_i = name_with_parenthesis(names{i});
    disp(name_i)
end
