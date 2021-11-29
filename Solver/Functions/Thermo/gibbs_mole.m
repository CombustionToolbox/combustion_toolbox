function val = gibbs_mole(mix)
    % Get the mole specific gibbs free energy [kJ/kmol]
    val = mix.g / mix.Ni * 1e3;
end