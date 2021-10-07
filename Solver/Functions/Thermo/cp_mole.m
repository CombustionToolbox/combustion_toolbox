function val = cp_mole(mix)
    % Get the mole-basis specific heat at constant pressure [kJ/kmol-K]
    val = mix.cP / mix.N;
end