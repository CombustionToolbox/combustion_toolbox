function val = cv_mass(mix)
    % Get the mass-basis specific heat at constant volume [kJ/kg-K]
    val = mix.cV / mix.mi * 1e-3;
end