function val = cv_mole(mix)
    % Get the mole-basis specific heat at constant volume [kJ/kmol-K]
    val = mix.cV / mix.N;
end