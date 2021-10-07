function val = cp_mass(mix)
    % Get the mass-basis specific heat at constant pressure [kJ/kg-K]
    val = mix.cP / mix.mi * 1e-3;
end