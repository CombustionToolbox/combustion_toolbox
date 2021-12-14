function g0 = set_g0(ls, TP, DB)
    % Function that computes the vector of gibbs free energy for the given
    % set of species [J/mol]
    for i=length(ls):-1:1
        species = ls{i};
        g0(i, 1) = species_g0(species, TP, DB) * 1e3;
    end
end