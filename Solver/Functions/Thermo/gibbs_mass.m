function val = gibbs_mass(mix)
    % Get the mass specific gibbs free energy [kJ/kg]
    val = mix.g / mix.mi;
end