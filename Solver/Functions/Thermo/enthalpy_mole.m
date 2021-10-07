function val = enthalpy_mole(mix)
    % Get the mole specific enthalpy [kJ/kmol]
    val = mix.h / mix.Ni * 1e3;
end